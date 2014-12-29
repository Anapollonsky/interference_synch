function [datout, pathout] = OFDMChan(data, channel)
siz = size(data);
datout = zeros(siz(1), siz(2));
pathout = zeros(siz(1), siz(2));

%% Random Configuration
    for k = 1:siz(2)
        [datout(:, k), pathout(:, k)] = step(channel, data(:, k));
%         reset(channel);
    end
end