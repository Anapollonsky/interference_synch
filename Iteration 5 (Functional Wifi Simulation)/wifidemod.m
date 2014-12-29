function [ datout ] = wifidemod(data, type, bpskdemod, qam4demod, qam16demod, qam64demod, ...
    condec16, condec78)
%wifidemod performs demodulation and decoding (potentially deinterleaving)
%on vector ata given certain rate. Assumes that certain objects exist.
    switch type
        case 1  %BPSK 1/2
            datademod = step(bpskdemod, data);
            datadec = step(condec16, real(datademod));
            
        case 2  %BPSK 3/4
            datademod = step(bpskdemod, data);
            datapar = reshape(datademod, 2, []);
            dataparfil = [datapar; zeros(1, length(data)/2)];
            datafil = reshape(dataparfil, [], 1);
            datadec = step(condec16, real(datafil));
            
        case 3  %QPSK 1/2
            datadec = step(condec16, data);
            datademod = step(qam4demod, data);
        case 4  %QPSK 3/4
            datademod = step(qam4demod, data);
            datadec = step(condec16, data);
            datadec(2:3:end) = [];
            
        case 5  %16QAM 1/2
            datademod = step(qam16demod, data);
            datadec = step(condec16, data);
            
           
        case 6  %16QAM 3/4
            datademod = step(qam16demod, data);
            datadec = step(condec16, data);
            datadec(2:3:end) = [];
            
        case 7  %64QAM 2/3
            datademod = step(qam64demod, data);
            datadec = step(condec78, data);
            datadec(3:4:end) = [];      
            
        case 8  %64QAM 3/4
            datademod = step(qam64demod, data);
            datadec = step(condec78, data);
            datadec(2:3:end) = [];        
            
    end
    datout = datadec;

% reset(errcalc16)
% data = randi([0 1], 3072, 1);
% dataenc = step(conenc16, data);
% % dataenc1 = dataenc;
% % dataenc1(3:3:end) = [];
% % datapar = reshape(dataenc1, 2, []);
% % dataparfil = [datapar; zeros(1, length(datamod)/2)];
% % datafil = reshape(dataparfil, [], 1);
% datafil = dataenc;
% datafil(4:6:end) = 0;
% datafil(5:6:end) = 0;
% datadec = step(condec16, real(datafil));
% step(errcalc16, data, datadec)
end
