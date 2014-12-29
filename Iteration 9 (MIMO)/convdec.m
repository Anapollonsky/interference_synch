function [datout] = convdec(data, condec, delay)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
    reset(condec);
    datadec = step(condec, data);
    datout = delayseq(datadec, -delay);
end

