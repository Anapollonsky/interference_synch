function [datout] = wifideint(data)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
global convdeint intdelay;
    release(convdeint); reset(convdeint);
    datadeint = step(convdeint, data);
    datout = delayseq(datadeint, -intdelay);
end

