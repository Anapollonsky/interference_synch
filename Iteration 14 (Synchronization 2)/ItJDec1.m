function [out1, out2, pg1, pg2] = ItJDec1(data, pg11, pg12, ...
    mod1, mod2, demod1, demod2)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
%     m1 = length(constellation(mod1));
    mu = .2;
    numiter = 50;
    out1 = zeros(length(data), 1);
    out2 = zeros(length(data), 1);
    pg1 = zeros(length(data), 1);
    pg1(1) = pg11;
    pg2 = zeros(length(data), 1); 
    pg2(1) = pg12;

    tempiter = 1;
    for t = 1:length(data)
%         sirl = pg1(tempiter) / pg2(tempiter);
        for k = 1:numiter
            %% Guess Symbol 1
            out1(t) = out1(t) + mu * ...
                ((data(t) - pg2(tempiter)*out2(t))/pg1(tempiter) - out1(t));
            %% Guess Symbol 2
            out2(t) = out2(t) + mu * ...
                ((data(t) - pg1(tempiter)*out1(t))/pg2(tempiter) - out2(t));
        end
        
        pg1(t) = pg1(tempiter) + mu * ...
            ((data(t) - pg2(tempiter)*out2(t))/out1(t) - pg1(tempiter));
        pg2(t) = pg2(tempiter) + mu * ...
            ((data(t) - pg1(tempiter)*out1(t))/out2(t) - pg2(tempiter));
        tempiter = t;
    end
end

