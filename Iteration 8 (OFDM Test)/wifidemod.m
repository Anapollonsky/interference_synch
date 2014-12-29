function [datout] = wifidemod(data, type)

    global bpskdemod qam4demod qam16demod qam64demod
%wifidemod performs demodulation based on transmission type.

    %% Assign modulation dimensions
    if type < 3
        demoder = bpskdemod;
    elseif type < 5
        demoder = qam4demod;
    elseif type < 7
        demoder = qam16demod;
    else
        demoder = qam64demod;
    end
%% Choose modulators based on dimensions

%     if m1 == 2
%         demoder = bpskdemod;
%     elseif m1 == 4
%         demoder = qam4demod;
%     elseif m1 == 16
%         demoder = qam16demod;
%     elseif m1 == 64
%         demoder = qam64demod;
%     end

%% Rest
    datout = step(demoder, data);
end