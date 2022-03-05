classdef BinaryStreamReader < handle
    % BINARYSTREAMREADER Memory-efficient streaming reader for raw radar frames.
    %
    % Reads binary radar captures in chunked frame iterations without loading
    % the entire continuous multi-gigabyte recording dataset into memory.
    
    properties
        FilePath                      % File path to binary recording
        FileId = -1;                  % Active file descriptor
        FrameSizeBytes                % Byte length of one frame
        TotalFrames                   % Total number of complete frames in file
        CurrentFrame = 0;             % 0-based frame counter
    end
    
    methods
        function obj = BinaryStreamReader(filePath, frameSizeBytes)
            % Construct a BinaryStreamReader instance.
            if nargin >= 1, obj.FilePath = filePath; end
            if nargin >= 2, obj.FrameSizeBytes = frameSizeBytes; end
            
            if ~isempty(obj.FilePath)
                obj.open();
            end
        end
        
        function delete(obj)
            obj.close();
        end
        
        function open(obj)
            % Open binary stream.
            obj.close();
            obj.FileId = fopen(obj.FilePath, 'rb');
            if obj.FileId < 0
                error('BinaryStreamReader:OpenFailed', 'Cannot open file: %s', obj.FilePath);
            end
            
            % Compute total frames
            fseek(obj.FileId, 0, 'eof');
            fileBytes = ftell(obj.FileId);
            fseek(obj.FileId, 0, 'bof');
            
            if isempty(obj.FrameSizeBytes) || obj.FrameSizeBytes <= 0
                obj.FrameSizeBytes = fileBytes;
            end
            
            obj.TotalFrames = floor(fileBytes / obj.FrameSizeBytes);
            obj.CurrentFrame = 0;
        end
        
        function close(obj)
            if obj.FileId >= 0
                fclose(obj.FileId);
                obj.FileId = -1;
            end
        end
        
        function [rawBytes, hasMore] = readNextFrame(obj)
            % Read next sequential frame bytes.
            if obj.FileId < 0 || obj.CurrentFrame >= obj.TotalFrames
                rawBytes = uint8([]);
                hasMore = false;
                return;
            end
            
            rawBytes = fread(obj.FileId, obj.FrameSizeBytes, '*uint8');
            obj.CurrentFrame = obj.CurrentFrame + 1;
            hasMore = (obj.CurrentFrame < obj.TotalFrames);
        end
        
        function seekFrame(obj, frameIndex)
            % Seek directly to 0-based frameIndex.
            if obj.FileId < 0
                error('BinaryStreamReader:NotOpen', 'Stream is not open.');
            end
            targetByte = frameIndex * obj.FrameSizeBytes;
            fseek(obj.FileId, targetByte, 'bof');
            obj.CurrentFrame = frameIndex;
        end
    end
end
