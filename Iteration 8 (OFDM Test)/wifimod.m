function [datout] = wifimod(data, type)
%wifimod Performs modulation based on tramission type.
global bpskmod qam4mod qam16mod qam64mod;
    if type < 3
        moder = bpskmod;
    elseif type < 5
        moder = qam4mod;
    elseif type < 7
        moder = qam16mod;
    else
        moder = qam64mod;
    end
    datout = step(moder, data);
end


