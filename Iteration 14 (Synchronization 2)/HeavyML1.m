function [out1, out2, error1, tapest2] = HeavyML1(data, d2, mod1, mod2, m1, m2, pg)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.
   
    tapest2 = zeros(length(data));
    error1 = zeros(length(data));
    frgt = .98;
    frgtvec = zeros(0);
    
    numsamp = 100;
    maxval = 1.5;
    tapsamp = [fliplr(-1 * logspace(-1, log10(maxval), numsamp/2)) ...
        logspace(-1, log10(maxval), numsamp/2)];
    
%     tapsamp = logspacespace(-maxval, maxval, numsamp);
    tapsind = 1:length(tapsamp);
    minimat = zeros(m1, numsamp, numsamp, length(data));
    maximat = zeros(m1, numsamp, numsamp, length(data));
    marximat = zeros(numsamp, numsamp, length(data));

    
    PG = squeeze(pg(1, :, :));
    
    for a = 1:m1-1
        aa = step(mod1, a);
        for r2 = tapsind
            r2s = tapsamp(r2);
            for i2 = tapsind
                i2s = tapsamp(i2);
                minimat(a+1, r2, i2, :) = ...
                    (PG) * aa;
            end
        end
    end
    
    for t = 1:length(data)
        bb = d2(t);
        for r2 = tapsind
            r2s = tapsamp(r2);
            for i2 = tapsind
                i2s = tapsamp(i2);
                minimat(:, r2, i2, t) =  minimat(:, r2, i2, t) + (r2s + 1i*i2s) * bb;
            end
        end

            
        maximat(:, :, :, t) = repmat(data(t), [m1 numsamp numsamp]) - minimat(:, :, :, t);
        
        for r2 = tapsind
            for i2 = tapsind
                [~, a] = min(maximat(:, r2, i2, t));
                marximat(r2, i2, t) = step(mod1, a-1);
            end
        end
        
        frgtvec = [frgtvec * frgt; 1];
        errcnt = mean( abs(  repmat(permute(d2(1:t), [2 3 1]), [numsamp numsamp 1]) - ...
                             marximat(:, :, 1:t) .* ...
                             repmat(permute(frgtvec, [2 3 1]), [numsamp numsamp 1])...
                          ), 3);
        errcnt2 = sum(sum(errcnt)) - errcnt;
        
        [err, minind] = max(reshape(errcnt2, [], 1));
        [a, b] = ind2sub(size(errcnt), minind);

        tapest2(t) = tapsamp(a) + 1i*tapsamp(b);
        error1(t) = err;
    end
    out1 = squeeze(marximat(a, b, :));
    out2 = d2; 
end

