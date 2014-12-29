function [datout] = wifimod(data, type)
%wifimod Performs encoding & modulation (potentially interleaving) on
%vector data given certain rate. Assumes that certain objects exist.
global convint bpskmod qam4mod qam16mod qam64mod conenc12 conenc23 conenc34 intdelay;
    switch type
        case 1  %BPSK 1/2
            release(convint);
            dataenc = step(conenc12, data);
            dataintp = step(convint, dataenc);
            dataint = delayseq(dataintp, intdelay);
            datamod = step(bpskmod, dataint);
        case 2  %BPSK 3/4
            release(convint);
            dataenc = step(conenc34, data);
            dataintp = step(convint, dataenc);
            dataint = delayseq(dataintp, intdelay);
            datamod = step(bpskmod, dataint);
        case 3  %QPSK 1/2
            release(convint);
            dataenc = step(conenc12, data);
            dataintp = step(convint, dataenc);
            dataint = delayseq(dataintp, intdelay);
            datamod = step(qam4mod, dataint);
        case 4  %QPSK 3/4
            release(convint);
            dataenc = step(conenc34, data);
            dataintp = step(convint, dataenc);
            dataint = delayseq(dataintp, intdelay);
            datamod = step(qam4mod, dataint);
        case 5  %16QAM 1/2
            release(convint);
            dataenc = step(conenc12, data);
            dataintp = step(convint, dataenc);
            dataint = delayseq(dataintp, intdelay);
            datamod = step(qam16mod, dataint);
        case 6  %16QAM 3/4
            release(convint);
            dataenc = step(conenc34, data);
            dataintp = step(convint, dataenc);
            dataint = delayseq(dataintp, intdelay);
            datamod = step(qam16mod, dataint);
        case 7  %64QAM 2/3
            release(convint);
            dataenc = step(conenc23, data);    
            dataintp = step(convint, dataenc);
            dataint = delayseq(dataintp, intdelay);
            datamod = step(qam64mod, dataint);
        case 8  %64QAM 3/4
            release(convint);
            dataenc = step(conenc34, data);
            dataintp = step(convint, dataenc);
            dataint = delayseq(dataintp, intdelay);
            datamod = step(qam64mod, dataint);
        case 9  %BPSK 1/1
            datamod = step(bpskmod, data);
    end
    datout = datamod;
end




