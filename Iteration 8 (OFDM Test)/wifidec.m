function [datout] = wifidec(data, type)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
    global condec12 condec23 condec34 delay12 delay23 delay34;
    reset(condec12);
    if type == 1 || type == 3 || type == 5
        condec = condec12;
        delay = delay12;
    elseif type == 7
        condec = condec23;
        delay = delay23;
    else
        condec = condec34;
        delay = delay34;
    end
    datadec = step(condec, data);
    datout = delayseq(datadec, -delay);
end

