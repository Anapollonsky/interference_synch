function [datout] = deinter(data, convint)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    release(convint); reset(convint);
    datout = step(convint, data);
end

