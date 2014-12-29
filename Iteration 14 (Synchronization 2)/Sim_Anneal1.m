function [out1, PG2, err] = Sim_Anneal1(data, d2, mod1, demod1, pg)
% SoftDecisionML1 Performs joint detection, including adaptive equalization.
iter = 10;
tempmax = 1000;
freqg = 0;
initg = 1;
assign = zeros(length(data), 1);


for k = 1:iter
    for temp = tempmax:-1:1
        
    end
end




end

function [pgout] = genpg(pg, 