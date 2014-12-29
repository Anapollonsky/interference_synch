function [datout] = wifienc(data, type)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
global conenc12 conenc23 conenc34;
    if type == 1 || type == 3 || type == 5
        conenc = conenc12;

    elseif type == 7
        conenc = conenc23;

    else
        conenc = conenc34;

    end
    datout = step(conenc, data);
end

