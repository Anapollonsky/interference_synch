function [out1, PG2, pgbounds] = IterativeML2(data, d2, mod1, pg)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
    m1 = length(constellation(mod1));
    numiter = 10;
    
    % This section controls forgetting factors. 
    keeplim = .1;
    frgt = .9;
    frgtkeep = ceil(log(keeplim)/log(frgt));
    frgtvec = zeros(length(data), length(data));
    frgkvec = zeros(numiter, numiter);
    frgk = .7;
    frgkkeep = ceil(log(keeplim)/log(frgk));
    frgkvec(1, 1) = 1;
    for k = 2:numiter, 
        frgkvec(1:k, k) = [frgkvec(1:k-1, k-1) * frgk; 1]; 
    end
    frgtvec(1, 1) = 1;
    for t = 2:length(data), 
        frgtvec(1:t, t) = [frgtvec(1:t-1, t-1) * frgt; 1]; 
    end
    
    out1 = zeros(length(data), 1);
       
    numsampPG = 50;
    pgmin = .1;
    PGSampMat = pgmatgen(-2, 2, -2, 2, pgmin, numsampPG);
    lastpgupdt = 1;
    pgbounds = repmat([-2 2 -2 2], [length(data) 1]); % xmin xmax ymin ymax
    
    minimat = zeros(numsampPG, numsampPG, length(data), numiter);
    maximat = zeros(m1, length(data), numiter);

    PG2 = zeros(length(data));
    PG = squeeze(pg(1));
    symguess = zeros(length(data));
    symconst = constellation(mod1);
    symconstep = .8;
   
    

    for t = 1:length(data)
%         PGSampMat
%         if t > 3         
%             if PG2(t-1) == PG2(t-2) && PG2(t-1) == PG2(t-3)
%                 pgbounds(t:end, :) = repmat([ ...
%                     mean([real(PG2(t-1)) pgbounds(t, 1)])...
%                     mean([real(PG2(t-1)) pgbounds(t, 2)]) ...
%                     mean([imag(PG2(t-1)) pgbounds(t, 3)]) ...
%                     mean([imag(PG2(t-1)) pgbounds(t, 4)]) ...
%                     ], [length(data) - t + 1 1]);     
%                 PGSampMat = pgmatgen(pgbounds(t, 1), pgbounds(t, 2), pgbounds(t, 3), ...
%                     pgbounds(t, 4), pgmin, numsampPG);
%     %             PGSampMat = PGSampMat;
%                 lastpgupdt = t;
%             end
%         end
        for k = 1:numiter
            %% Guess at PG2, assume sent modulated 0
            minimat(:, :, t, k) = (data(t) - PG * symguess(t) - PGSampMat * d2(t)).^2;

            miniavg = mean( mean( abs( minimat(:, :, max(lastpgupdt, t-frgtkeep+1):t, max(1, k-frgkkeep+1):k) .*  ...
                repmat(permute(frgtvec(max(lastpgupdt, t-frgtkeep+1):t, t), [2 3 1]), [numsampPG, numsampPG, 1, min(k, frgkkeep)]) .* ...
                repmat(permute(frgkvec(max(1, k-frgkkeep+1):k, k), [2 3 4 1]), [numsampPG, numsampPG, min(t - lastpgupdt + 1, frgtkeep), 1])...
                ), 3), 4);
            temp = reshape(miniavg, [], 1);
            [~, minind] = min(temp);
            [a, b] = ind2sub(size(miniavg), minind);
            
            %% Guess at modulated, using calculated PG2
            maximat(:, t, k) = (data(t) - PG * symconst - PGSampMat(a, b) * d2(t)).^2;
            maxiavg = mean( abs( maximat(:, t, max(1, k-frgkkeep+1):k) ...
                .* repmat(permute(frgkvec(max(1, k-frgkkeep+1):k, k), [2 3 1]), [m1, 1, 1]) ...
                ), 3);
            [~, maxind] = min(maxiavg);
            symguess(t) = symguess(t) + ...
                (symconst(maxind) - symguess(t)) * symconstep;
        end
        out1(t) = symguess(t);
        PG2(t) = PGSampMat(a, b); 
    end
end

function [out] = pgmatgen(xmin, xmax, ymin, ymax, minmag, num)
    logxmin = min([max(abs(xmax), minmag) max(abs(xmin), minmag)]);
    logxmax = max(abs(xmax), abs(xmin));
    logymin = min([max(abs(ymax), minmag) max(abs(ymin), minmag)]);
    logymax = max(abs(ymax), abs(ymin));
    if abs(xmax + xmin) < logxmax
        infactx = log10(logxmin) - log10(minmag);
        outfactx = log10(logxmax) - log10(logxmin);
        inumx = floor((num * (2 * infactx / (2 * infactx + outfactx))) / 2);
        outnumx = num - 2 * inumx;
        if abs(xmax) > abs(xmin)
            tapsampPGX = [fliplr(...
                -1 * logspace(log10(minmag), log10(logxmin), inumx)) ...
                logspace(log10(minmag), log10(logxmax), inumx + outnumx)];
        else
            tapsampPGX = [fliplr(...
                -1 * logspace(log10(minmag), log10(logxmax), outnumx + inumx))...
                logspace(log10(minmag), log10(logxmin), inumx)];  
        end
    else
        if sign(xmax ) == 1
            tapsampPGX = logspace(max(log10(minmag), log10(xmin)), log10(xmax), num);
        else
            tapsampPGX = fliplr(...
                -1 * logspace(max(log10(minmag), log10(abs(xmax))), log10(abs(xmin)), num));
        end
    end
    
    if abs(ymax + ymin) < logymax
        infacty = log10(logymin) - log10(minmag);
        outfacty = log10(logymax) - log10(logymin);
        inumy = floor((num * (2 * infacty / (2 * infacty + outfacty))) / 2);
        outnumy = num - 2 * inumy;
        
        if abs(xmax) > abs(xmin)
            tapsampPGY = [fliplr(...
                -1 * logspace(log10(minmag), log10(logymin), inumy)) ...
                logspace(log10(minmag), log10(logymax), inumy + outnumy)];
        else
            tapsampPGY = [fliplr(...
                -1 * logspace(log10(minmag), log10(logymax), outnumy + inumy))...
                logspace(log10(minmag), log10(logymin), inumy)];  
        end
    else
        if sign(ymax ) == 1
            tapsampPGY = logspace(max(log10(minmag), log10(ymin)), log10(ymax), num);
        else
            tapsampPGY = fliplr(...
                -1 * logspace(max(log10(minmag), log10(abs(ymax))), log10(abs(ymin)), num));
        end
    end
    
    [X, Y] = meshgrid(1:num, 1:num);
    out = tapsampPGX(X) + 1i* tapsampPGY(Y);
end

