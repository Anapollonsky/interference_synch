function [perr] = sicserrqfunc(snr, inr, sir)
% sicserrqfunc calculates the symbol error for the successive inteference cancellation estimator for
% a single signal with the Q function.
    if sir > 4
        perr = 1/2 * qfunc(sqrt(snr) - 2 * sqrt(inr));
    elseif sir > 9/4
        perr = 1/2 * (1-qfunc(sqrt(snr) - sqrt(inr)));
    elseif sir > 1
        perr = 1/2 * (1 - qfunc(2*sqrt(inr) - sqrt(snr)));
    elseif sir > 1/4
        perr = 1/2 * qfunc(sqrt(inr) - sqrt(snr));
    else
        perr = qfunc(sqrt(snr));
    end
end
