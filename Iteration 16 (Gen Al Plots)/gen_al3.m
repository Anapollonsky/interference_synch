function [out1, PG2, iters, pgbase, rotfac] = gen_al3(data, d2, mod1, pg)
% Genetic algorithm approach to decoding primary signal while synchronizing
% interferer.
const = constellation(mod1);
m = length(const);
len = length(data);

% Convergence Checking
symbthresh = 1;
analogvarthresh = 1;

% Other parameters
numsols = 60; % 4x4 -> 40
numcheck = round(sqrt(numsols));
maxpg = 3;
minpg = .02;
maxfreq = pi/10; % Randomly picks from entire circle of possible phase change per symbol.

% Initialization
iter = true;
iters = 1;
sols = zeros(numsols, len + 2);
rotvec = zeros(numsols^2 + numsols, len);
combos = zeros(len+2, numsols, numsols);
for k = 1:numsols % Initial PG, Phase-Change-Per-Turn, symbols
    sols(k, :) = [newpg(minpg, maxpg) newang(maxfreq) const(randi(m, [1 len])).'];
end

% Genetic Stuff
while iter == true
    % Cross-mix the current solution vector with itself
    for p = 1:numsols
        for o = 1:numsols
            ranum = rand(len + 2, 1);
            for k = 1:2
                if ranum(k) < .25
                    combos(k, p, o) = sols(p, k);
                elseif ranum(k) < .5
                    combos(k, p, o) = sols(o, k);
                elseif ranum(k) < .7
                    combos(k, p, o) = mean([sols(p, k) sols(o, k)]);
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
                    combos(k, p, o) = sols(p, k);
                elseif ranum(k) < .9
                    combos(k, p, o) = sols(o, k);
                else
                    combos(k, p, o) = const(randi(m, 1));
                end
            end
        end
    end
    comboser = [reshape(combos, 2 + len, [], 1) sols.']; % Combos reshaped into linear form + old values
    for t = 1:len
        rotvec(:, t) = exp(1i*comboser(2, :)).^(t-1);
    end
    
    % Evaluating fitness and picking best resulting 2-D solution matrix
    fits = mean((abs(repmat(data, 1, numsols^2 + numsols) - ...
        comboser(3:end, :) .* repmat(pg, 1, numsols^2 + numsols) - ...
        repmat(comboser(1, :), len, 1) .* rotvec.' .* ...
        repmat(d2, 1, numsols^2 + numsols))), 1).';
   
    [~, bestind] = sort(fits, 1, 'ascend');
    sols = comboser(:, bestind(1:numsols)).';
    
    % Checking for convergence
    best = sols(1:numcheck, :);
    decmod = repmat(mode(best(:, 3:end), 1), numcheck, 1);
    decindx = mean(mean(decmod == best(:, 3:end)));
    angvar = var(best(:, 2));
    pgvar = var(best(:, 1));
    if decindx >= symbthresh && angvar <= analogvarthresh && pgvar <= analogvarthresh
        iter = false;
    else
        iters = iters + 1;
    end
end
out1 = sols(1, 3:end).';
PG2 = rotvec(1, :).* sols(1, 1);
rotfac = sols(1, 2);
pgbase = sols(1, 1);
end

function [pgout] = newpg(ming, maxg)
    pgout = (ming + (maxg-ming).*rand(1, 1)) * exp(1i * 2*pi*rand(1, 1));
end

function [angout] = newang(maxrot)
    angout = -maxrot + (2*maxrot * rand(1, 1));
end

            
        