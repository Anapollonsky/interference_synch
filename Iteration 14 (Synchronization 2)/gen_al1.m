function [out1, PG2] = gen_al1(data, d2, mod1, demod1, pg)
% SoftDecisionML1 Performs joint detection, including adaptive equalization.
len = length(data);
iter = 5;
iter2 = 5;
numsol = 250;
newpick = 80;
maxpg = 3;
minpg = .1;
maxfreq = pi/5;

assigns = zeros(len, iter);
assignsmod = zeros(len, iter);
sols1 = zeros(numsol, 2);
rotvec = zeros(numsol, len);
rotvec2 = zeros(newpick^2 + newpick, len);
rotvec3 = zeros(len, 1);
combos = zeros(2, newpick, newpick);
for k = 1:numsol % Initial PG Magnitude, PG Angle, Phase-Change-Per-Turn
    sols1(k, :) = [(minpg + (maxpg-minpg).*rand(1, 1)) ...
        * exp(1i * 2*pi*rand(1, 1)) ...
        -maxfreq + (2*maxfreq * rand(1, 1))];
end

for k = 1:iter
    for l = 1:iter2
        for t = 1:len
            rotvec(:, t) = exp(1i*sols1(:, 2)).^(t-1);
        end

        fits = mean(abs(repmat(data, 1, numsol) - repmat(assignsmod(:, k), 1, numsol)...
            .* repmat(pg, 1, numsol) - ...
            repmat( permute(sols1(:, 1), [2 1]), len, 1) ...
            .* rotvec.' .* repmat(d2, 1, numsol)), 1).'; % numsol x 1 distance

        [~, bestind] = sort(fits, 1, 'ascend');
        best = sols1(bestind(1:newpick), :);
        for p = 1:newpick
            for o = 1:newpick
                ranum = rand(2, 1);
                if ranum(1) < .25
                    combos(1, p, o) = best(p, 1);
                elseif ranum(1) < .5
                    combos(1, p, o) = best(o, 1);
                elseif ranum(1) < .85
                    combos(1, p, o) = mean([best(p, 1) best(o, 1)]);
                else
                    combos(1, p, o) = (minpg + (maxpg-minpg).*rand(1, 1)) ...
        * exp(1i * 2*pi*rand(1, 1));
                end

                if ranum(2) < .25
                    combos(2, p, o) = best(p, 2);
                elseif ranum(2) < .5
                    combos(2, p, o) = best(o, 2);
                elseif ranum(2) < .85
                    combos(2, p, o) = mean([best(p, 2) best(o, 2)]);
                else
                    combos(2, p, o) = -maxfreq + (2*maxfreq * rand(1, 1));
                end
            end
        end
        combos2 = [reshape(combos, 2, [], 1) best.'];

        for t = 1:len
            rotvec2(:, t) = exp(1i*combos2(2, :)).^(t-1);
        end

        fits2 = mean(abs(repmat(data, 1, newpick^2 + newpick) - ...
            repmat(assignsmod(:, k), 1, newpick^2+newpick) .* repmat(pg, 1, newpick^2+newpick) - ...
            repmat(combos2(1, :), len, 1) ...
            .* rotvec2.' .* repmat(d2, 1, newpick^2+newpick)), 1).'; % numsol x 1 distance

        [~, bestind2] = sort(fits2, 1, 'ascend');
        sols1 = combos2(:, bestind2(1:numsol)).';
    end
    bestsol = sols1(1, :).';
    for t = 1:len
        rotvec3(t) = exp(1i*bestsol(2)).^(t-1);
    end
    assigns(:, k) = (data - rotvec3.* sols1(1) .* d2) ./ pg;
    assignsmod(:, k) = step(mod1, step(demod1, assigns(:, k)));
end

out1 = assignsmod(:, end);
PG2 = rotvec3.* sols1(1);

end
