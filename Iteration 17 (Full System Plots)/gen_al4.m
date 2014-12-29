function [out1, PG2, iters, pgbase, rotfac] = gen_al4(data, d2, mod1, pg)
% Genetic algorithm approach to decoding primary signal while synchronizing
% interferer.
const = constellation(mod1);
m = length(const);
len = length(data);

% Convergence Checking
symbthresh = 1;
analogvarthresh = .001;

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
% combos = zeros(len+2, numsols, numsols);
for k = 1:numsols % Initial PG, Phase-Change-Per-Turn, symbols
    sols(k, :) = [newpg(minpg, maxpg, 1, 1) newang(maxfreq, 1, 1) const(randi(m, [1 len])).'];
end

% sols = [newpg(minpg, maxpg, numsols, 1) newang(maxfreq, numsols, 1) const(randi(m, [len numsols ]))]
% Genetic Stuff
while iter == true
    % Cross-mix the current solution vector with itself
    ranum = rand(len+2, numsols, numsols);
    combos = zeros(len+2, numsols, numsols);
    sols2 = permute(sols, [3 1 2]); % 1 numsols len+2
    for k = 1:2
        ranum2 = ranum(k, :, :) < .25;
        combos(k, :, :) = combos(k, :, :) + ranum2 .* repmat(sols2(:, :, k), [1 1 numsols]);
        
        ranum2 = .25 < ranum(k, :, :) & ranum(k, :, :) < .5;
        combos(k, :, :) = combos(k, :, :) + ranum2 .* repmat(permute(sols2(:, :, k), [1 3 2]), [1 numsols 1]);
        
        ranum2 = .5 < ranum(k, :, :) & ranum(k, :, :) < .7;
        combos(k, :, :) = combos(k, :, :) ...
            + ranum2 .*(repmat(permute(sols2(:, :, k), [1 3 2]), [1 numsols 1]) ...
            + repmat(sols2(:, :, k), [1 1 numsols])) / 2;
        
        ranum2 = .7 < ranum(k, :, :);
        if k == 1
            combos(k, :, :) = combos(k, :, :) + ranum2 .* permute(newpg(minpg, maxpg, numsols, numsols), [3 1 2]);
        else
            combos(k, :, :) = combos(k, :, :) + ranum2 .* permute(newang(maxfreq, numsols, numsols), [3 1 2]);
        end
    end
    for k = 3:len+2
        ranum2 = ranum(k, :, :) < .45;
        combos(k, :, :) = combos(k, :, :) + ranum2 .* repmat(sols2(:, :, k), [1 1 numsols]);
        ranum2 = .45 < ranum(k, :, :) & ranum(k, :, :) < .9;
        combos(k, :, :) = combos(k, :, :) + ranum2 .* repmat(permute(sols2(:, :, k), [1 3 2]), [1 numsols 1]);
        ranum2 = .9 < ranum(k, :, :);
        combos(k, :, :) = combos(k, :, :) + ranum2 .* permute(const(randi(m, numsols)), [3 1 2]);              
    end     
%% old    
%     for p = 1:numsols
%         for o = 1:numsols
%             ranum = rand(len + 2, 1);
%             for k = 1:2
%                 if ranum(k) < .25
%                     combos(k, p, o) = sols(p, k);
%                 elseif ranum(k) < .5
%                     combos(k, p, o) = sols(o, k);
%                 elseif ranum(k) < .7
%                     combos(k, p, o) = mean([sols(p, k) sols(o, k)]);
%                 else
%                     if k == 1
%                         combos(k, p, o) = newpg(minpg, maxpg);
%                     else
%                         combos(k, p, o) = newang(maxfreq);
%                     end
%                 end
%             end
%             for k = 3:len+2
%                 if ranum(k) < .45
%                     combos(k, p, o) = sols(p, k);
%                 elseif ranum(k) < .9
%                     combos(k, p, o) = sols(o, k);
%                 else
%                     combos(k, p, o) = const(randi(m, 1));
%                 end
%             end
%         end
%     end
%% ctd    
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
    
    if iters == 100
        break;
    end
end
out1 = sols(1, 3:end).';
PG2 = rotvec(bestind(1), :).* sols(1, 1);
rotfac = sols(1, 2);
pgbase = sols(1, 1);
end

function [pgout] = newpg(ming, maxg, y, x)
    pgout = (ming + (maxg-ming).*rand(y, x)) * exp(1i * 2*pi*rand(y, x));
end

function [angout] = newang(maxrot, y, x)
    angout = -maxrot + (2*maxrot * rand(y, x));
end

            
        