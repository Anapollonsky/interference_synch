function [out1, PG2] = gen_al2(data, d2, mod1, pg)
% SoftDecisionML1 Performs joint detection, including adaptive equalization.
const = constellation(mod1);
m = length(const);
len = length(data);
iter = 100; % ~10 for 4/4
numsol = 5000;
newpick = 200;
maxpg = 3;
minpg = .02;
maxfreq = pi;

sols1 = zeros(numsol, len + 2);
rotvec = zeros(numsol, len);
rotvec2 = zeros(newpick^2 + newpick, len);
rotvec3 = zeros(len, 1);
combos = zeros(len+2, newpick, newpick);
for k = 1:numsol % Initial PG, Phase-Change-Per-Turn, symbols
    sols1(k, :) = [newpg(minpg, maxpg) newang(maxfreq) const(randi(m, [1 len])).'];
end


for l = 1:iter
    for t = 1:len
        rotvec(:, t) = exp(1i*sols1(:, 2)).^(t-1);
    end

    fits = mean(abs(repmat(data, 1, numsol) - sols1(:, 3:end).' .* ...
        repmat(pg, 1, numsol) - ...
        repmat( permute(sols1(:, 1), [2 1]), len, 1) ...
        .* rotvec.' .* repmat(d2, 1, numsol)), 1).';

    [~, bestind] = sort(fits, 1, 'ascend');
    best = sols1(bestind(1:newpick), :);
    for p = 1:newpick
        for o = 1:newpick
            ranum = rand(len + 2, 1);
            for k = 1:2
                if ranum(k) < .3
                    combos(k, p, o) = best(p, k);
                elseif ranum(k) < .6
                    combos(k, p, o) = best(o, k);
                elseif ranum(k) < .9
                    combos(k, p, o) = mean([best(p, k) best(o, k)]);
                else
                    if k == 1
                        combos(k, p, o) = newpg(minpg, maxpg);
                    else
                        combos(k, p, o) = newang(maxfreq);
                    end
                end
            end
            for k = 3:len+2
                if ranum(k) < .45
                    combos(k, p, o) = best(p, k);
                elseif ranum(k) < .9
                    combos(k, p, o) = best(o, k);
                else
                    combos(k, p, o) = const(randi(m, 1));
                end
            end
        end
    end
    combos2 = [reshape(combos, 2 + len, [], 1) best.'];
    for t = 1:len
        rotvec2(:, t) = exp(1i*combos2(2, :)).^(t-1);
    end
    
    
    fits2 = mean((abs(repmat(data, 1, newpick^2 + newpick) - ...
        combos2(3:end, :) .* repmat(pg, 1, newpick^2 + newpick) - ...
        repmat(combos2(1, :), len, 1) .* rotvec2.' .* ...
        repmat(d2, 1, newpick^2 + newpick))), 1).';
   
    [~, bestind2] = sort(fits2, 1, 'ascend');
    sols1 = combos2(:, bestind2(1:numsol)).';
end
bestsol = sols1(1, :).';
for t = 1:len
    rotvec3(t) = exp(1i*bestsol(2)).^(t-1);
end
out1 = bestsol(3:end);
PG2 = rotvec3.* sols1(1);
end


function [pgout] = newpg(ming, maxg)
    pgout = (ming + (maxg-ming).*rand(1, 1)) * exp(1i * 2*pi*rand(1, 1));
end

function [angout] = newang(maxrot)
    angout = -maxrot + (2*maxrot * rand(1, 1));
end

            
        