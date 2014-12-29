function [datout] = qrm_mld(data, PG, mod, M)
% Implementation of QRM_MLD. Only works with 2x2 MIMO.
    const = constellation(mod);
    datout = zeros(length(data), 2);
    for t = 1:length(data) % Loop over time
        [Q, R] = qr(squeeze(PG(t, :, :)).'); % QR Decomposition
        yhat = Q'*(data(t, :).'); % Generate Q*y
        dist = zeros(length(const), 1); % Distance Preallocation
        for a = 1:length(const)
            dist(a) = abs(yhat(2) - R(2, 2)*const(a)); % Calculate distance
        end
        [~, I] = sort(dist, 'ascend'); %Find best solutions
        x2s = const(I(1:min(length(const), M))); % Grab best solutions, move to stage 2
        dist = zeros(length(x2s), length(const)); % Distance Preallocation
        for a = 1:length(x2s) % Loop over last iteration
            for b = 1:length(const) % Loop over current possibilities
                dist(a, b) = (yhat(2) - R(2, 2)*x2s(a))^2 + ...
                    (yhat(1) - R(1, 1)*const(b) - R(1, 2)*x2s(a))^2; % Calculate distance
            end
        end
        [~, minind] = min(reshape(dist, [], 1)); % Find best index
        [a, b] = ind2sub(size(dist), minind); % Find corresponding symbols
        datout(t, :) = [const(b) x2s(a)]; % Output them
    end
end

