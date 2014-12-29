classdef GoldEnc < matlab.System
%MIMOEnc separates data and spreads it across several antennas. Time-rows
%are interleaved. Assumes properly sized data frame.
    methods (Access=protected)
        function out = stepImpl(~, data)
            theta = (1+sqrt(5))/2;
            thetahat = 1-theta;
            alp = 1 + 1i*(1-theta);
            alphat = 1+1i*(1-thetahat);
            out = zeros(length(data)/2, 2);
            for k = 1:length(data)/4
                datatemp = data((k-1)*4 + 1: 4*k);
                out((k-1)*2 + 1:2*k, :) = 1/sqrt(5) * [...
                    alp*(datatemp(1) + theta*datatemp(2)) ...
                    alp *(datatemp(3) + theta*datatemp(4));...
                    1i * alphat * (datatemp(3) + thetahat * datatemp(4))...
                    alphat * (datatemp(1) + thetahat*datatemp(2))...
                    ];
            end
        end
        
        function numIn = getNumInputsImpl(~)
           numIn = 1;
        end
 
        function numOut = getNumOutputsImpl(~)
           numOut = 1;
        end
    end
end