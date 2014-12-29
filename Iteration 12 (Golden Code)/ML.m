classdef ML < matlab.System
%Uses ML for decoding FDFR.
    methods (Access=protected)
        function out = stepImpl(~, data, PG)
            global bpskmod
            theta = (1+sqrt(5))/2;
            thetahat = 1-theta;
            alp = 1 + 1i*(1-theta);
            alphat = 1 + 1i*(1-thetahat);
            m = 2;
            maximat = zeros(m, m, m, m, length(data)/2);
            minimat1 = zeros(length(data), 2, 2, m, m, m, m);
            minimat2 = zeros(length(data), 2, 2, m, m, m, m);
            out = zeros(length(data)*2, 1);
            for t = 1:2:length(data)
                tempg1 = [PG(t, 1, 1) PG(t, 2, 1);...
                          PG(t+1, 1, 1) PG(t+1, 2, 1)];
                tempg2 = [PG(t, 1, 2) PG(t, 2, 2);...
                          PG(t+1, 1, 2) PG(t+1, 2, 2)];
                for a = 0:m-1
                    aa = step(bpskmod, a);
                    for b = 0:m-1
                        bb = step(bpskmod, b);
                        for c = 0:m-1
                            cc = step(bpskmod, c);
                            for d = 0:m-1
                                dd = step(bpskmod, d);
                                tempcode = [alp *(aa + theta * bb) alp * (cc + theta*dd); ...
                                    1i*alphat * (cc + thetahat*dd) alphat * (aa + thetahat*bb)];
                                
                                minimat1(t, :, :, a+1, b+1, c+1, d+1) = permute(...
                                    1/sqrt(5) * tempg1 .* tempcode, [3 1 2]);
                                
                                minimat2(t, :, :, a+1, b+1, c+1, d+1) = permute(...
                                    1/sqrt(5) * tempg2 .* tempcode, [3 1 2]);
                                
                                maximat(a+1, b+1, c+1, d+1, (t+1)/2) = ...
                                    sum(sum(abs((data(t:t+1, 1) - sum(permute(minimat1(t, :, :, a+1, b+1, c+1, d+1), [2 3 1]), 2)) ) )) + ...
                                    sum(sum(abs((data(t:t+1, 2) - sum(permute(minimat2(t, :, :, a+1, b+1, c+1, d+1), [2 3 1]), 2)) ) ));
                            end
                        end
                    end
                end
                [~, minind] = min(reshape(maximat(:, :, :, :, (t+1)/2), [], 1));
                [aaa, bbb, ccc, ddd] = ind2sub(size(maximat(:, :, :, :, 1)), minind); 
                out(((t-1)*2 + 1):(2*(t+1))) = [aaa-1; bbb-1; ccc-1; ddd-1];    
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