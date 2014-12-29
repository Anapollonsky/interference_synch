classdef ZFMIMODec < matlab.System
% Andrew Apollonsky. Implementation of Zero-Forcing linear MIMO decoder. 
    methods (Access=protected)
        function out = stepImpl(~, data, PG, numTx)
            out = zeros(length(data), numTx);
            W = zeros(numTx, numTx, length(data));
            for t = 1:length(data)
                k = squeeze(PG(t, :, :));
                W(:, :, t) = (k'*k)\k';
                out(t, :) = data(t, :)*W(:, :, t);
            end
        end
        
        function numIn = getNumInputsImpl(~)
           numIn = 3;
        end
 
        function numOut = getNumOutputsImpl(~)
           numOut = 1;
        end
    end
end