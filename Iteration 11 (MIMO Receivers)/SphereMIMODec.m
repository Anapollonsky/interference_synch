classdef JMDMIMODec < matlab.System
%JMDMIMODEC uses JMD, an ML-like detector, to decode MIMO signals. Accepts
%up to 8x8 MIMO, does not accept non-square configurations. 
    methods (Access=protected)
        function out = stepImpl(~, data, PG, mod, m, numTx, numRx)
            % Matrix preallocation
            out = zeros(length(data), numTx);
            mvec = [ones(numRx, 1)*m; ones(8-numRx, 1)];
            minimat1 = zeros(length(data), numTx, numRx, m);
            minimat2 = zeros(length(data), numRx, mvec(1), mvec(2), mvec(3), mvec(4), mvec(5), mvec(6), mvec(7), mvec(8));
            grid1 = ones(length(data), mvec(1), mvec(2), mvec(3), mvec(4), mvec(5), mvec(6), mvec(7), mvec(8), numTx);             
            
            % Generation of values for JMD detector
            for p = 1:numTx
                for n = 0:m-1
                    aa = step(mod, n);
                    for k = 1:numRx
                        minimat1(:, p, k, n+1) = PG(:, p, k) * aa;
                    end
                end
            end
            
            % Grid for repmatting later.
            gridthing(1, :) = [1 2 3 4 5 6 7 8 9];
            for k = 2:8
                gridthing(k, :) = gridthing(k-1, :);
                gridthing(k, [k+1, k]) = gridthing(k, [k, k+1]); 
            end
            
            % Generating values that will be subtracted for JMD receiver.
            for k = 1:numRx
                for p = 1:numTx
                    grid1(:, :, :, :, :, :, :, :, :, p) = repmat(permute(squeeze(minimat1(:, p, k, :)), gridthing(p, :)), [1 ones(1, p-1)*m 1 ones(1, numTx-p)*m]);
                end
                minimat2(:, k, :, :, :, :, :, :, :, :) = permute(sum(grid1, 10), [1 10 2 3 4 5 6 7 8 9]);
            end
          
           
            % Calculating Maximat, which stores the differences between the
            % actual data and would-be data for every single combination.
            maximat = permute(squeeze(sum(abs(squeeze(repmat(data, [1, 1, mvec(1), mvec(2), mvec(3), mvec(4), mvec(5), mvec(6), mvec(7), mvec(8)])) - ...
                                      squeeze(minimat2)), 2)), [2 3 4 5 6 7 8 9 1]);
                                  
            % Using the values to find the decoded signal.
            for t = 1:length(data)
                [~, minind] = min(reshape(maximat( :, :, :, :, :, :, :, :, t), [], 1));
                siz = size(maximat(:, :, :, :, :, :, :, :, t));
                [a, b, c, d, e, f, g, h] = ind2sub(siz, minind);
                k = [a b c d e f g h];
                out(t, :) = k(1:numTx)-1;               
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
