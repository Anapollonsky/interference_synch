function [out1, PG2] = IterativeML1(data, d2, mod1, pg)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
    m1 = length(constellation(mod1));
    numiter = 10;
    
    frgt = .95;
    frgtvec = zeros(length(data), length(data));
    frgkvec = zeros(numiter, numiter);
    frgk = .7;
    frgkvec(1, 1) = 1;
    for k = 2:numiter, 
        frgkvec(1:k, k) = [frgkvec(1:k-1, k-1) * frgk; 1]; 
    end
    
    frgtvec(1, 1) = 1;
    for t = 2:length(data), 
        frgtvec(1:t, t) = [frgtvec(1:t-1, t-1) * frgt; 1]; 
    end
    
    out1 = zeros(length(data), 1);
       
    numsampPG = 10;
    maxvalPG = 2;
    tapsampPG = [fliplr(-1 * logspace(-1, log10(maxvalPG), numsampPG/2)) ...
        logspace(-1, log10(maxvalPG), numsampPG/2)];
    [X, Y] = meshgrid(1:numsampPG, 1:numsampPG);
    PGSampMat = tapsampPG(X) + 1i* tapsampPG(Y);
   
    minimat = zeros(numsampPG, numsampPG, length(data), numiter);
    maximat = zeros(m1, length(data), numiter);

    PG2 = zeros(length(data));
    PG = squeeze(pg(1));
    symguess = zeros(length(data));
    symconst = constellation(mod1);
    symconstep = .3;
   

    for t = 1:length(data)      
        for k = 1:numiter
            %% Guess at PG2, assume sent modulated 0
            minimat(:, :, t, k) = (data(t) - PG * symguess(t) - PGSampMat * d2(t)).^2;
            miniavg = mean( mean( abs( minimat(:, :, 1:t, 1:k) .*  ...
                repmat(permute(frgtvec(1:t, t), [2 3 1]), [numsampPG, numsampPG, 1, k]) .* ...
                repmat(permute(frgkvec(1:k, k), [2 3 4 1]), [numsampPG, numsampPG, t, 1])...
                ), 3), 4);
            temp = reshape(miniavg, [], 1);
            [~, minind] = min(temp);
            [a, b] = ind2sub(size(miniavg), minind);
            
            %% Guess at modulated, using calculated PG2
            maximat(:, t, k) = (data(t) - PG * symconst - PGSampMat(a, b) * d2(t)).^2;
            maxiavg = mean( abs( maximat(:, t, 1:k) ...
                .* repmat(permute(frgkvec(1:k, k), [2 3 1]), [m1, 1, 1]) ...
                ), 3);
            [~, maxind] = min(maxiavg);
            symguess(t) = symguess(t) + ...
                (symconst(maxind) - symguess(t)) * symconstep;
        end
        out1(t) = symguess(t);
        PG2(t) = PGSampMat(a, b); 
    end
end

