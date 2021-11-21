classdef Crc32
    % CRC32 IEEE 802.3 standard 32-bit Cyclic Redundancy Check calculator.
    %
    % Used for validating integrity of TI DCA1000 EVM raw LVDS packets and
    % serialized binary radar cube frames.
    
    properties (Constant, Access = private)
        % Standard IEEE 802.3 reversed polynomial
        POLYNOMIAL = uint32(hex2dec('EDB88320'));
    end
    
    methods (Static)
        function crc = compute(dataBytes)
            % Compute 32-bit CRC checksum over an array of bytes or characters.
            %
            % Input:
            %   dataBytes - uint8 array or string of ASCII bytes
            % Output:
            %   crc       - uint32 CRC value
            
            if ischar(dataBytes)
                dataBytes = uint8(dataBytes);
            else
                dataBytes = uint8(dataBytes(:));
            end
            
            table = chirp.utils.Crc32.getTable();
            crc = uint32(hex2dec('FFFFFFFF'));
            
            poly = chirp.utils.Crc32.POLYNOMIAL;
            len = length(dataBytes);
            
            for k = 1:len
                byteVal = uint32(dataBytes(k));
                idx = bitxor(bitand(crc, uint32(255)), byteVal) + 1;
                crc = bitxor(bitshift(crc, -8), table(idx));
            end
            
            crc = bitxor(crc, uint32(hex2dec('FFFFFFFF')));
        end
        
        function hexStr = computeHex(dataBytes)
            % Compute 32-bit CRC and format as an 8-character uppercase hex string.
            val = chirp.utils.Crc32.compute(dataBytes);
            hexStr = sprintf('%08X', val);
        end
    end
    
    methods (Static, Access = private)
        function table = getTable()
            % Generate or retrieve static 256-element lookup table.
            persistent crcTable
            if isempty(crcTable)
                crcTable = zeros(256, 1, 'uint32');
                poly = uint32(hex2dec('EDB88320'));
                for byteVal = 0:255
                    c = uint32(byteVal);
                    for bit = 1:8
                        if bitand(c, uint32(1))
                            c = bitxor(bitshift(c, -1), poly);
                        else
                            c = bitshift(c, -1);
                        end
                    end
                    crcTable(byteVal + 1) = c;
                end
            end
            table = crcTable;
        end
    end
end
