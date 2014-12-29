function [out1, PG2, iters] = gen_al1(data, d2, mod1, pg)
% Genetic algorithm approach to decoding primary signal while synchronizing
% interferer.
const = constellation(mod1);
m = length(const);
len = size(data, 1);
numofdm = size(data, 2);

% Convergence Checking
symbthresh = 1;
analogvarthresh = 1;

% Other parameters
numsols = 100; % 4x4 -> 40
numcheck = round(sqrt(numsols));
maxpg = 3;
minpg = .02;
maxfreq = pi; % Randomly picks from entire circle of possible phase change per symbol.

% Initialization
iter = true;
iters = 1;
sols = zeros(numsols, numofdm * (len + 1) + 1);
rotvec = zeros(len, numofdm, numsols^2 + numsols);
combos = zeros(numofdm * (len + 1) + 1, numsols, numsols);
for k = 1:numsols % Initial PG, Phase-Change-Per-Turn, symbols
    sols(k, 1) = newang(maxfreq);
    for p = 1:numofdm
        sols(k, p+1) = newpg(minpg, maxpg);
    end
    sols(k, numofdm + 2 : end) = const(randi(m, [1 len*numofdm])).';
end

% Genetic Stuff
while iter == true
    % Cross-mix the current solution vector with itself
    for p = 1:numsols
        for o = 1:numsols
            ranum = rand(numofdm * (len + 1) + 1, 1);
            % Frequency offset
            if ranum(1) < .25
                combos(1, p, o) = sols(p, 1);
            elseif ranum(1) < .5
                combos(1, p, o) = sols(o, 1);
            elseif ranum(1) < .7
                combos(1, p, o) = mean([sols(p, 1) sols(o, 1)]);
            else
                combos(1, p, o) = newang(maxfreq);
            end
            % Path Gains
            for k = 2:numofdm + 1
                if ranum(k) < .25
                    combos(k, p, o) = sols(p, k);
                elseif ranum(k) < .5
                    combos(k, p, o) = sols(o, k);
                elseif ranum(k) < .7
                    combos(k, p, o) = mean([sols(p, k) sols(o, k)]);
                else
                    combos(k, p, o) = newpg(minpg, maxpg);
                end
            end
            % Data Estimates
            for k = numofdm + 2:numofdm * (len + 1) + 1 
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
    comboser = [reshape(combos, numofdm * (len + 1) + 1 , [], 1) sols.']; % Combos reshaped into linear form + old values
    
    combodec = reshape(comboser(numofdm + 2:end, :), len, numofdm, []);
    combopg = reshape(comboser(2:numofdm + 1, :), [1 numofdm numsols^2 + numsols]);
    comboang = reshape(comboser(1, :), 1, 1, []);
    for t = 1:len
        rotvec(t, :, :) = repmat(exp(1i*comboang).^(t-1), [1 numofdm]);
    end

    % Evaluating fitness and picking best resulting 2-D solution matrix
    fits = squeeze( mean( mean( abs( repmat(data, [1 1 numsols^2 + numsols]) - ...
        combodec .* repmat(pg, [1 1 numsols^2 + numsols]) - ...
        repmat(combopg, [len 1 1]) .* ...
        rotvec .* ...
        repmat(d2, [1 1 numsols^2 + numsols])), 1 ), 2));
   
    [~, bestind] = sort(fits, 1, 'ascend');
    sols = comboser(:, bestind(1:numsols)).';
    
    % Checking for convergence
    
    best = sols(1:numcheck, :);
    decmod = repmat(mode(best(:, 2 + numofdm:end), 1), numcheck, 1);
    decindx = mean(mean(decmod == best(:, 2 + numofdm:end)));
    angvar = var(best(:, 1));
    pgvar = mean(var(best(:, 2:numofdm+1), 1), 2);
    if decindx >= symbthresh && angvar <= analogvarthresh && pgvar <= analogvarthresh
        iter = false;
    else
        iters = iters + 1;
    end
end
out1 = reshape(sols(1, numofdm + 2:end), len, numofdm);
PG2 = rotvec(:, :, bestind(1)).* repmat(sols(1, 2:numofdm + 1), len, 1);
end

function [pgout] = newpg(ming, maxg)
    pgout = (ming + (maxg-ming).*rand(1, 1)) * exp(1i * 2*pi*rand(1, 1));
end

function [angout] = newang(maxrot)
    angout = -maxrot + (2*maxrot * rand(1, 1));
end

            
        