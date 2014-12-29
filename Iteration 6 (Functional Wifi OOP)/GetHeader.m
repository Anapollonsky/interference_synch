function [TransType ModOrder Rate BinLength] = GetHeader(data)
%UNTITLED7 Takes the 13th time point of a frame and derives signal
%characteristics.
    demodata = wifidemod(data.', 1).';
    ttbin = demodata(1:3);
    k = mat2str(ttbin);
    ttbinstr = strcat(k(2), k(4), k(6));
    TransType = bin2dec(ttbinstr) + 1;

    if TransType < 3
        ModOrder = 2;
    elseif TransType < 5
        ModOrder = 4;
    elseif TransType < 7
        ModOrder = 16;
    else
        ModOrder = 64;
    end

    if TransType == 1 || TransType == 3 || TransType == 5
        Rate = 1/2;
    elseif TransType == 7
        Rate = 2/3;
    else
        Rate = 3/4;
    end

    lenbin = demodata(4:19);
    p = mat2str(lenbin);
    lenbinstr = strcat(
    
end

