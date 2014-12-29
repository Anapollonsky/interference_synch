classdef MMSESICMIMODec < matlab.System
%MMSESICMIMODEC uses MMSE-SIC to decode MIMO signals. Supports 8x8 MIMO. 
    methods (Access=protected)
        function out = stepImpl(~, data, PG, N0, mod, demod, numTx)
            out = zeros(length(data), numTx);  
            pow = sum(PG.^2, 3);
            [~, powind] = sort(pow, 2, 'descend');
            for t = 1:length(data)
                PGT = squeeze(PG(t, :, :)).'; 
                data2 = data(t, :);
                for k = 1:numTx
                    W = (PGT'*PGT + eye(numTx)*N0)\PGT';
                    temp = W*data2.';
                    xhat = temp(powind(t, k));
                    out(t, powind(t, k)) = step(demod, xhat);
                    mods = step(mod, out(t, powind(t, k)));
                    data2 = data2 - PGT(:, powind(t, k)).' * mods;
                end
            end
        end
        
        function numIn = getNumInputsImpl(~)
           numIn = 6;
        end      
        
        function numOut = getNumOutputsImpl(~)
           numOut = 1;
        end
    end
end