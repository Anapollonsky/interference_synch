function [datout] = convenc(data, conenc)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
    datout = step(conenc, data);
end

