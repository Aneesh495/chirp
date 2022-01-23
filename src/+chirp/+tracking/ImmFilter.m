classdef ImmFilter < handle
    % IMMFILTER Interacting Multiple Model estimator for maneuvering radar targets.
    %
    % Dynamically switches between multiple sub-filter models (e.g. constant velocity,
    % aggressive acceleration) via a discrete Markov mode transition matrix.
    
    properties
        Models                        % Cell array of sub-filters (e.g. {ekf_cv, ekf_ca})
        TransitionMatrix              % N_models x N_models Markov transition probability matrix Pi
        ModelProbabilities            % 1 x N_models mode probability vector mu
    end
    
    methods
        function obj = ImmFilter(filterModels, transitionMatrix, initialProbabilities)
            % Construct an ImmFilter instance.
            if nargin >= 1, obj.Models = filterModels; end
            if nargin >= 2 && ~isempty(transitionMatrix)
                obj.TransitionMatrix = transitionMatrix;
            else
                % Default 2-model transition matrix
                obj.TransitionMatrix = [0.95, 0.05; 0.10, 0.90];
            end
            if nargin >= 3 && ~isempty(initialProbabilities)
                obj.ModelProbabilities = initialProbabilities(:).';
            else
                nMod = length(obj.Models);
                obj.ModelProbabilities = ones(1, nMod) / nMod;
            end
        end
        
        function [mixedStates, mixedCovs] = computeMixingProbabilities(obj, states, covs)
            % Step 1: Model interaction / state mixing.
            nMod = length(obj.Models);
            Pi = obj.TransitionMatrix;
            mu = obj.ModelProbabilities;
            
            % Normalizing constants c_j = sum_i (Pi_ij * mu_i)
            c = mu * Pi;
            
            % Mixing probabilities: mu_ij = (Pi_ij * mu_i) / c_j
            mu_mix = zeros(nMod, nMod);
            for i = 1:nMod
                for j = 1:nMod
                    mu_mix(i, j) = (Pi(i, j) * mu(i)) / max(c(j), 1e-12);
                end
            end
            
            dim = length(states{1});
            mixedStates = cell(1, nMod);
            mixedCovs = cell(1, nMod);
            
            for j = 1:nMod
                % Mixed state mean x0_j = sum_i (mu_ij * x_i)
                x0_j = zeros(dim, 1);
                for i = 1:nMod
                    x0_j = x0_j + mu_mix(i, j) * states{i};
                end
                mixedStates{j} = x0_j;
                
                % Mixed covariance P0_j = sum_i (mu_ij * [P_i + (x_i - x0_j)(x_i - x0_j)'])
                P0_j = zeros(dim, dim);
                for i = 1:nMod
                    diffX = states{i} - x0_j;
                    P0_j = P0_j + mu_mix(i, j) * (covs{i} + diffX * diffX');
                end
                mixedCovs{j} = P0_j;
            end
        end
        
        function [combinedState, combinedCov, updatedStates, updatedCovs] = step(...
                obj, states, covs, dt, measurement)
            % Execute a complete IMM cycle: mixing, prediction, update, and mode probability update.
            
            nMod = length(obj.Models);
            [mixedStates, mixedCovs] = obj.computeMixingProbabilities(states, covs);
            
            updatedStates = cell(1, nMod);
            updatedCovs = cell(1, nMod);
            likelihoods = zeros(1, nMod);
            
            for j = 1:nMod
                filter = obj.Models{j};
                [pState, pCov] = filter.predict(mixedStates{j}, mixedCovs{j}, dt);
                [uState, uCov, innov, S] = filter.update(pState, pCov, measurement);
                
                updatedStates{j} = uState;
                updatedCovs{j} = uCov;
                
                % Gaussian likelihood: Lambda_j = N(innov; 0, S)
                mDim = length(innov);
                detS = max(det(S), 1e-12);
                invS = inv(S);
                quadForm = innov' * invS * innov;
                normConst = 1.0 / sqrt(((2.0 * pi)^mDim) * detS);
                likelihoods(j) = max(normConst * exp(-0.5 * quadForm), 1e-15);
            end
            
            % Step 3: Model probability update
            % mu_j = c_j * Lambda_j / sum(c_k * Lambda_k)
            c = obj.ModelProbabilities * obj.TransitionMatrix;
            unnormMu = c .* likelihoods;
            obj.ModelProbabilities = unnormMu / max(sum(unnormMu), 1e-12);
            
            % Step 4: Output state combination
            dim = length(updatedStates{1});
            combinedState = zeros(dim, 1);
            for j = 1:nMod
                combinedState = combinedState + obj.ModelProbabilities(j) * updatedStates{j};
            end
            
            combinedCov = zeros(dim, dim);
            for j = 1:nMod
                diffX = updatedStates{j} - combinedState;
                combinedCov = combinedCov + obj.ModelProbabilities(j) * ...
                    (updatedCovs{j} + diffX * diffX');
            end
        end
    end
end
