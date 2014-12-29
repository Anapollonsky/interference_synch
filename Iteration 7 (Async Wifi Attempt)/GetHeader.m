function [TransType, ModOrder, Rate, BinLength] = GetHeader(bindata)
%GetHeader Takes the 13-15th time points of a frame and derives signal
%characteristics.
    % Demodulate and all that
%     demodata = wifidemod(reshape(data, [], 1), 1);
%     bindata = wifidec(wifideint(demodata), 1);
    
    % Derive Transmission Type to find TransType, ModOrder and Rate
    ttbin = bindata(1:3);
    k = mat2str(ttbin);
    ttbinstr = strcat(k(2), k(4), k(6));
    TransType = bin2dec(ttbinstr) + 1;

    %Mod Order
    if TransType < 3
        ModOrder = 2;
    elseif TransType < 5
        ModOrder = 4;
    elseif TransType < 7
        ModOrder = 16;
    else
        ModOrder = 64;
    end
    
    % Code Rate
    if TransType == 1 || TransType == 3 || TransType == 5
        Rate = 1/2;
    elseif TransType == 7
        Rate = 2/3;
    else
        Rate = 3/4;
    end

    % Find Binary Length
    lenbin = bindata(4:19);
    p = mat2str(lenbin);
    lenbinaccu = '';
    for iter = 1:16
        indx = 2:2:32;
        lenbinaccu = strcat(lenbinaccu, p(indx(iter)));
    end
    BinLength = bin2dec(lenbinaccu);
end

