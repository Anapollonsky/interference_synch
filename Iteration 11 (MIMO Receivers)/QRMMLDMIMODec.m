classdef QRMMLDMIMODec < matlab.System
%JMDMIMODEC uses JMD, an ML-like detector, to decode MIMO signals. Accepts
%up to 8x8 MIMO, does not accept non-square configurations. 

    methods (Access=protected)
        function out = stepImpl(~, data, PG, numTx, mod, M2)
            const = constellation(mod);
%             bitable = [0;1];
            out = zeros(length(data), numTx);
%             if length(const) < M2
%                 M = length(const);
%             else
%                 M = M2;
%             end
            M = M2;
            accudat = zeros(length(const), M);
            choic = zeros(numTx, M);
            choosemat1 = zeros(M);
            choosemat2 = zeros(M);
            for t = 1:length(data)
                [Q, R] = qr(squeeze(PG(t, :, :)).');
                yhat = Q'*(data(t, :).');
%                 temp = (yhat(numTx) - R(numTx, numTx)*const).^2;
%                 [~, best] = sort(temp, 1, 'ascend');
% %                 dist(1, 1:M) = temp(best(1:M)).';
%                 choic(1, 1:M) = temp(best(1:M)); 
                
                for k = 1:numTx
                    for o = 1:length(const)
                        for p = 1:M % Generate error matrix
                            accudat(o, p) = sum(...
                                (yhat(numTx - k + 1:end) ...
                                - R(numTx - k + 1:end, numTx - k + 1:end) ...
                                * [const(o); choic(1:k-1, p)]).^2 ...
                                , 1);
                        end
                    end
                    if k ~= numTx  % If not the last iteration, use sort
                        [~, best] = sort(reshape(accudat, [], 1), 'ascend');
                        choic2 = choic;
                        for p = 1:M
                            [choosemat1(p), choosemat2(p)]...
                                = ind2sub(size(accudat), best(p));
                            choic(:, p) = [const(choosemat1(p));...
                                choic2(1:end-1, choosemat2(p))];
                        end
                    else % If last iteration, pick best value
                        [~, best] = min(reshape(accudat, [], 1));
                        [choosemat1(1), choosemat2(1)] = ind2sub(size(accudat), best);
                        out(t, :) = [const(choosemat1(1))...
                                choic(1:end-1, choosemat2(1)).'];
                    end
                end

            end
        end
        
        function numIn = getNumInputsImpl(~)
           numIn = 5;
        end
 
        function numOut = getNumOutputsImpl(~)
           numOut = 1;
        end
    end
end
