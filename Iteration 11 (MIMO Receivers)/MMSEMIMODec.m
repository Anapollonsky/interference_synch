classdef MMSEMIMODec < matlab.System
%MIMOEnc separates data and spreads it across several antennas. Time-rows
%are interleaved. Assumes properly sized data frame.
    methods (Access=protected)
        function out = stepImpl(~, data, PG, N0, numTx)
            out = zeros(length(data), numTx);
            W = zeros(numTx, numTx, length(data));
            for t = 1:length(data)
                k = squeeze(PG(t, :, :));
                W(:, :, t) = (k'*k + eye(numTx)*N0)\k';
                out(t, :) = data(t, :)*W(:, :, t);
            end
        end
        
        function numIn = getNumInputsImpl(~)
           numIn = 4;
        end
 
        function numOut = getNumOutputsImpl(~)
           numOut = 1;
        end
    end
end