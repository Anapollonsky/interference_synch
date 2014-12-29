function [perr] = intigserr(m1,sinr)
% jointmdserr calculates the symbol error for the interference ignorant
% detector.
perr = 2 * (1 - 1/(sqrt(m1))) * (m1 - 1) ./ (3*sinr);
end

