function [out1, PG2, iters] = gen_al2(data, d2, mod1, pg)
% Genetic Algorithm approach to decoding OFDM symbols while synchronizing
% a same-subcarrier interferer.

const = constellation(mod1);
m = length(const);
len = size(data, 1);
numofdm = size(data, 2);

% Convergence Check
symbthresh = 1;
analogvarthresh = .0005;

% Other parameters
numsols = 100;
numcheck = round(sqrt(numsols));
maxpg = 3;
minpg = .02;
maxfreq = pi;

% Initialization
iter = true;
iters = 1;
combos = zeros(numofdm * (len + 1) + 1, numsols, numsols);
sols = zeros(numsols, numofdm * (len + 1) + 1);
rotvec = zeros(len, numofdm, numsols^2 + numsols);

for k = 1:numsols
    sols(k, 1) = newang(maxfreq);
    for p = 1:numofdm
        sols(k, p+1) = newpg(minpg, maxpg);
    end
    sols(k, numofdm + 2 :end) = const(randi(m, [1 len*numofdm])).';
end

while iter == true
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
    comboser = [reshape(combos, numofdm * (len + 1) + 1, [], 1,
            







end

function [pgout] = newpg(ming, maxg)
    pgout = (ming + (maxg-ming).*rand(1, 1)) * exp(1i * 2*pi*rand(1, 1));
end

function [angout] = newang(maxrot)
    angout = -maxrot + (2*maxrot * rand(1, 1));
end
