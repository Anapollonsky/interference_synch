classdef OFDMPost < matlab.System
%OFDMPost Process Time x Subcarrier signal for OFDM transmission.
    methods (Access=protected)
        function out = stepImpl(~, data, mu)
            fd = data(mu+1:end, :);
            out = fftshift(fft(fd)).';
        end
    end
end