function [datout] = inter(data, convdeint, intdelay)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
    release(convdeint); reset(convdeint);
    datadeint = step(convdeint, data);
    datout = delayseq(datadeint, -intdelay);
end

