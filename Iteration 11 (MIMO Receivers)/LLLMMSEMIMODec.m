classdef LLLMMSEMIMODec < matlab.System
%MMSESICMIMODEC uses MMSE-SIC to decode MIMO signals. Supports 8x8 MIMO. 
    methods (Access=protected)
        function out = stepImpl(~, data, PG, N0)
            del = .25;
            T = eye(2);
            out = zeros(length(data), 2);
            for t = 1:length(data)
                    [Q, R] = qr(squeeze(PG(t, :, :)).');
%                     yhat = Q'*data(t, :)

                for numiter = 1:100     
                    
                    mu12 = round(R(1, 2)/R(1, 1));
                    R(:, 2) = R(:, 2) - mu12 * R(:, 1);
                    T(:, 2) = T(:, 2) - mu12 * T(:, 1);
                    
                    if del * R(1, 1)^2 > R(1, 2)^2 + R(2, 2)^2
                        R = [R(:, 2) R(:, 1)];
                        T = [T(:, 2) T(:, 1)];
                        a1 = R(1, 1) / sqrt(R(1, 1)^2 + R(2, 1)^2);
                        b1 = R(2, 1) / sqrt(R(1, 1)^2 + R(2, 1)^2);
                        phi1 = [a1 b1; -b1 a1];
                        R = phi1*R;
                        Q = Q*phi1.';
                        continue;
                    end
                    %% hide
%                     mu12 = round(R(1, 2)/R(1, 1));
%                     R(:, 2) = R(:, 2) - mu12 * R(:, 1);
%                     T(:, 2) = T(:, 2) - mu12 * T(:, 1);
%                     
%                     if del * R(1, 1)^2 > R(1, 2)^2 + R(2, 2)^2
%                         R = [R(:, 2) R(:, 1) R(:, 3:4)];
%                         T = [T(:, 2) T(:, 1) T(:, 3:4)];
%                         a1 = R(1, 1) / sqrt(R(1, 1)^2 + R(2, 1)^2);
%                         b1 = R(2, 1) / sqrt(R(1, 1)^2 + R(2, 1)^2);
%                         phi1 = [a1 b1 0 0; -b1 a1 0 0; 0 0 1 0; 0 0 0 1];
%                         R = phi1*R;
%                         Q = Q*phi1.';
%                         continue;
%                     end



%                     mu23 = round(R(2, 3) / R(2, 2));
%                     R(:, 3) = R(:, 3) - mu23*R(:, 2);
%                     T(:, 3) = T(:, 3) - mu23*T(:, 2);
%                     mu13 = round(R(1, 3) / R(1, 1));
%                     R(:, 3) = R(:, 3) - mu13*R(:, 1);
%                     T(:, 3) = T(:, 3) - mu13*T(:, 1);
%                     
%                     if del*R(2, 2) > R(2, 3)^2 + R(3, 3)^2
%                         R = [R(:, 1) R(:, 3) R(:, 2) R(:, 4)];
%                         T = [T(:, 1) T(:, 3) T(:, 2) T(:, 4)];
%                         a2 = R(2, 2) / sqrt(R(2, 2)^2 + R(3, 2)^2);
%                         b2 = R(3, 2) / sqrt(R(2, 2)^2 + R(3, 2)^2);
%                         phi2 = [1 0 0 0; 0 a2 b2 0; 0 -b2 a2 0; 0 0 0 1];
%                         R = phi2 * R;
%                         Q = Q*phi2.';
%                         continue;
%                     end
%                     
%                     mu34 = round(R(3, 4)/R(3, 3));
%                     R(:, 4) = R(:, 4) - mu34*R(:, 3);
%                     T(:, 4) = T(:, 4) - mu34*T(:, 3);
%                     mu24 = round(R(2, 4)/R(2, 2));
%                     R(:, 4) = R(:, 4) - mu24*R(:, 2);
%                     T(:, 4) = T(:, 4) - mu24*T(:, 2);
%                     mu14 = round(R(1, 4)/R(1, 1));
%                     R(:, 4) = R(:, 4) - mu14*R(:, 1);
%                     T(:, 4) = T(:, 4) - mu14*T(:, 1);
%                     
%                     if del*R(3, 3) > R(3, 4)^2 + R(4, 4)^2
%                         R = [R(:, 1:2) R(:, 4) R(:, 3)];
%                         T = [T(:, 1:2) T(:, 4) T(:, 3)];
%                         a3 = R(3, 3) / sqrt(R(3, 3)^2 + R(4, 3)^2);
%                         b3 = R(4, 3) / sqrt(R(3, 3)^2 + R(4, 3)^2);
%                         phi3 = [1 1 0 0; 1 1 0 0; 0 0 a3 b3; 0 0 -b3 a3];
%                         R = phi3 * R;
%                         Q = Q*phi3.';
%                         continue;
%                     end
%% 
                    break;
                end
                yhat = Q'*data(t, :).';
                xhat = (R'*R + eye(2)*N0)\R'*yhat;
                xhat = round(xhat);
%                 size(T)
%                 size(xhat)
                out(t, :) = (T\xhat).';
                
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