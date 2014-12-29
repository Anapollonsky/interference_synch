classdef OFDMPre < matlab.System
%OFDMPre Pre-Process Time x Subcarrier signal for OFDM transmission.
    methods (Access=protected)
        function out = stepImpl(~, data)
            FFTPoints = length(data');
            td = ifft(data, FFTPoints, 1);
            out = [td(1:16, :); td(:, :)];
        end
    end
end