classdef MIMOEnc < matlab.System
%MIMOEnc separates data and spreads it across several antennas. Time-rows
%are interleaved. Assumes properly sized data frame.
    methods (Access=protected)
        function out = stepImpl(~, data, numTx)
            siz = size(data);
            out = zeros(siz(1)/numTx, siz(2), numTx);
            for k = 1:numTx
                out(:, :, k) = data(k:numTx:end, :);
            end
        end
        
        function numIn = getNumInputsImpl(~)
           numIn = 2;
        end
 
        function numOut = getNumOutputsImpl(~)
           numOut = 1;
        end
    end
end