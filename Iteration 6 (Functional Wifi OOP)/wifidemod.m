function [datout] = wifidemod(data, type)

    global bpskdemod qam4demod qam16demod qam64demod condec12 ...
            condec23 condec34 convdeint
%wifidemod performs demodulation and decoding (potentially deinterleaving)
%on vector ata given certain rate. Assumes that certain objects exist.

    %% Assign modulation dimensions
    if type < 3
        m1 = 2;
    elseif type < 5
        m1 = 4;
    elseif type < 7
        m1 = 16;
    else
        m1 = 64;
    end
    
    %% Choose modulators based on dimensions

    if m1 == 2
        demoder = bpskdemod;
    elseif m1 == 4
        demoder = qam4demod;
    elseif m1 == 16
        demoder = qam16demod;
    elseif m1 == 64
        demoder = qam64demod;
    end

    %% Choose Decoder
    if type == 1 || type == 3 || type == 5
        condec = condec12;
    elseif type == 7
        condec = condec23;
    else
        condec = condec34;
    end
    
%     reset(convdeint);
    release(convdeint);
    
%     release(condec);
    datademod = step(demoder, data);
    datadeint = step(convdeint, datademod);
%     datadeint = datademod;
    datout = step(condec, datadeint);
%     reset(condec);
   
end