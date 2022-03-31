function test_no_em_dashes()
    % TEST_NO_EM_DASHES Verification audit ensuring no em dashes exist anywhere in repo.
    %
    % Scans all source code, documentation, scripts, and build manifests to assert
    % that the Unicode em dash character (U+2014) is never used.
    
    fprintf('Running test_no_em_dashes...\n');
    
    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(thisDir));
    
    % Recursively find all tracked text files
    files = getAllFiles(projectRoot);
    
    % Unicode em dash sequence in UTF-8
    emDashChar = char([226, 128, 148]);
    
    violationCount = 0;
    
    for k = 1:length(files)
        filePath = files{k};
        
        % Skip .git, binary output, and build directories
        if ~isempty(strfind(filePath, '.git')) || ~isempty(strfind(filePath, 'bin/mex')) || ...
           ~isempty(strfind(filePath, 'output')) || ~isempty(strfind(filePath, '.rcube')) || ...
           ~isempty(strfind(filePath, '.mex')) || ~isempty(strfind(filePath, '.o'))
            continue;
        end
        
        try
            content = fileread(filePath);
        catch
            continue;
        end
        
        matches = strfind(content, emDashChar);
        if ~isempty(matches)
            fprintf('  VIOLATION: Em dash found in file: %s (count: %d)\n', filePath, length(matches));
            violationCount = violationCount + length(matches);
        end
    end
    
    fprintf('  Total em dash occurrences found: %d\n', violationCount);
    assert(violationCount == 0, 'Em dash policy violation: Em dashes are strictly prohibited.');
    fprintf('  PASSED: test_no_em_dashes (Zero em dashes verified across codebase)\n');
end

function fileList = getAllFiles(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    
    if ~isempty(fileList)
        fileList = cellfun(@(f) fullfile(dirName, f), fileList, 'UniformOutput', false);
    end
    
    subDirs = {dirData(dirIndex).name};
    validIndex = ~ismember(subDirs, {'.', '..', '.git'});
    
    for iDir = find(validIndex)
        nextDir = fullfile(dirName, subDirs{iDir});
        fileList = [fileList; getAllFiles(nextDir)];
    end
end
