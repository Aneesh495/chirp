classdef RadarCubeCodec
    % RADARCUBECODEC High-throughput binary serialization for multi-channel radar cubes.
    %
    % Writes and reads structured .rcube binary files with magic byte validation,
    % dimension headers, float64 payloads, and CRC32 integrity checksums.
    
    properties (Constant)
        MAGIC_HEADER = 'CHIRP_CUBE_V1';
    end
    
    methods (Static)
        function save(filePath, radarCube, metadata)
            % Serialize 3D radar cube to .rcube binary file with CRC32 checksum.
            
            if nargin < 3, metadata = struct(); end
            
            fid = fopen(filePath, 'wb');
            if fid < 0
                error('RadarCubeCodec:FileOpenFailed', 'Cannot create file: %s', filePath);
            end
            cleanUp = onCleanup(@() fclose(fid));
            
            % 1. Write Magic Header (13 bytes)
            fwrite(fid, chirp.io.RadarCubeCodec.MAGIC_HEADER, 'char');
            
            % 2. Write Dimensions [N_adc, N_chirps, N_virt] as uint32
            dims = uint32(size(radarCube));
            if length(dims) < 3, dims(3) = 1; end
            fwrite(fid, dims(1:3), 'uint32');
            
            % 3. Write Complex Payload (single or double precision)
            % Interleave real and imaginary components as float64
            nElements = prod(dims(1:3));
            payload = zeros(2 * nElements, 1, 'double');
            payload(1:2:end) = real(radarCube(:));
            payload(2:2:end) = imag(radarCube(:));
            
            fwrite(fid, payload, 'double');
            
            % 4. Compute and append CRC32 checksum
            byteData = typecast(payload, 'uint8');
            crcVal = chirp.utils.Crc32.compute(byteData);
            fwrite(fid, crcVal, 'uint32');
        end
        
        function [radarCube, dims] = load(filePath)
            % Deserialize and validate a .rcube binary file.
            
            fid = fopen(filePath, 'rb');
            if fid < 0
                error('RadarCubeCodec:FileNotFound', 'Cannot open file: %s', filePath);
            end
            cleanUp = onCleanup(@() fclose(fid));
            
            % 1. Verify Magic Header
            headerLen = length(chirp.io.RadarCubeCodec.MAGIC_HEADER);
            header = fread(fid, headerLen, '*char').';
            if ~strcmp(header, chirp.io.RadarCubeCodec.MAGIC_HEADER)
                error('RadarCubeCodec:InvalidMagicHeader', 'File does not contain valid CHIRP_CUBE_V1 header.');
            end
            
            % 2. Read Dimensions
            dims = fread(fid, 3, 'uint32').';
            
            % 3. Read Payload
            nElements = prod(dims);
            payload = fread(fid, 2 * nElements, 'double');
            
            % 4. Read CRC32 trailer
            storedCrc = fread(fid, 1, 'uint32');
            
            % Verify CRC32
            byteData = typecast(payload, 'uint8');
            calcCrc = chirp.utils.Crc32.compute(byteData);
            if ~isempty(storedCrc) && (storedCrc ~= calcCrc)
                error('RadarCubeCodec:ChecksumMismatch', ...
                    'CRC32 checksum verification failed! Stored: %08X, Computed: %08X', ...
                    storedCrc, calcCrc);
            end
            
            % Reconstruct complex array
            re = payload(1:2:end);
            im = payload(2:2:end);
            complexLinear = re + 1j * im;
            
            radarCube = reshape(complexLinear, dims);
        end
    end
end
