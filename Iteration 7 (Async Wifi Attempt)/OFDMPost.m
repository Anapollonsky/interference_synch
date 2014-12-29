classdef OFDMPost < matlab.System
%OFDMPost Process Time x Subcarrier signal for OFDM transmission.
    methods (Access=protected)
        function out = stepImpl(~, data)
            fd = data(16+1:end, :);
            FFTPoints = length(fd');
            out = fft(fd, FFTPoints, 1);
        end
    end
end