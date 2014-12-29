function [out1, out2, error1, tapest1, tapest2] = HeavyML1(data, d1, d2, mod1, mod2, m1, m2)
%SoftDecisionML1 Performs joint detection, including adaptive equalization.

    k = size(zeros(m1, m2));
    
    tapest1 = zeros(length(data));
    tapest2 = zeros(length(data));
    error1 = zeros(length(data));
    
    numsamp = 20;
    maxval = 1.5;
    tapsamp = [fliplr(-1 * logspace(-1, log10(maxval), numsamp/2)) ...
        logspace(-1, log10(maxval), numsamp/2)];
    
%     tapsamp = logspacespace(-maxval, maxval, numsamp);
    tapsind = 1:length(tapsamp);
    minimat = zeros(m1, m2, numsamp, numsamp, numsamp, numsamp);
    maximat = zeros(m1, m2, numsamp, numsamp, numsamp, numsamp, length(data));
    marximat = zeros(numsamp, numsamp, numsamp, numsamp, length(data), 2);

    for r1 = tapsind
        r1s = tapsamp(r1);
        for i1 = tapsind
            i1s = tapsamp(i1);
            for r2 = tapsind
                r2s = tapsamp(r2);
                for i2 = tapsind
                    i2s = tapsamp(i2);
                    for a = 0:m1-1
                        aa = step(mod1, a);
                        for b = 0:m2-1
                            bb = step(mod2, b);
                            minimat(a+1, b+1, r1, i1, r2, i2) = ...
                                (r1s + 1i*i1s) * aa + (r2s + 1i*i2s) * bb;
                        end
                    end
                end
            end
        end
    end
    
    for t = 1:length(data)
        maximat(:, :, :, :, :, :, t) = repmat(data(t), [m1 m2 numsamp numsamp numsamp numsamp]) - minimat;
        
        for r1 = tapsind
            for i1 = tapsind
                for r2 = tapsind
                    for i2 = tapsind
                        [~, minind] = min(reshape(maximat(:, :, r1, i1, r2, i2, t), [], 1));
                        [a, b] = ind2sub(k, minind);
                        marximat(r1, i1, r2, i2, t, 1) = step(mod1, a-1);
                        marximat(r1, i1, r2, i2, t, 2) = step(mod2, b-1);
                    end
                end
            end
        end   
        
        errcnt = mean( (repmat(permute(d1(1:t), [2 3 4 5 1]), [numsamp numsamp numsamp numsamp 1]) - marximat(:, :, :, :, 1:t, 1)) .^2 ...
                     + (repmat(permute(d2(1:t), [2 3 4 5 1]), [numsamp numsamp numsamp numsamp 1]) - marximat(:, :, :, :, 1:t, 2)) .^2 , 5);
           
        [err, minind] = min(reshape(errcnt, [], 1));
        [a, b, c, d] = ind2sub(size(errcnt), minind);

        tapest1(t) = tapsamp(a) + 1i*tapsamp(b);
        tapest2(t) = tapsamp(c) + 1i*tapsamp(d);
        error1(t) = err;
    end
    out1 = squeeze(marximat(a, b, c, d, :, 1));
    out2 = squeeze(marximat(a, b, c, d, :, 2)); 
end

