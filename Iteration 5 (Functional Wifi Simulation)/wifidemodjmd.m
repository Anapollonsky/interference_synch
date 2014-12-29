function [datout1, datout2] = wifidemodjmd(data, type1, type2, g1, g2, ...
    bpskmod, qam4mod, qam16mod, qam64mod,...
    bpskdecmod, qam4decmod, qam16decmod, qam64decmod,...
    bpskdemod, qam4demod, qam16demod, qam64demod,...
    condec12, condec23, condec34, convdeint)
%wifidemod performs demodulation and decoding (potentially deinterleaving)
%on vector ata given certain rate. Assumes that certain objects exist.

    %% Assign modulation dimensions
    if type1 < 3
        m1 = 2;
    elseif type1 < 5
        m1 = 4;
    elseif type1 < 7
        m1 = 16;
    else
        m1 = 64;
    end
    if type2 < 3
        m2 = 2;
    elseif type2 < 5
        m2 = 4;
    elseif type2 < 7
        m2 = 16;
    else
        m2 = 64;
    end     
    
    %% Choose modulators based on dimensions
    % There are unnecessary modulations due to the JMD algorithm using a
    % non- bitinput/bitoutput modulator, while the rest of the script
    % working in bits. Thus, once decoded by JMD, the data is reincoded
    % with the normal modulator and decoded with the bitoutput demodulator.
    if m1 == 2
        umod1 = bpskdecmod;
        remod1 = bpskmod;
        demod1 = bpskdemod;
    elseif m1 == 4
        umod1 = qam4decmod;
        remod1 = qam4decmod;
        demod1 = qam4demod;
    elseif m1 == 16
        umod1 = qam16decmod;
        remod1 = qam16decmod;
        demod1 = qam16demod;
    elseif m1 == 64
        umod1 = qam64decmod;
        remod1 = qam64decmod;
        demod1 = qam64demod;
    end
    
    if m2 == 2
        umod2 = bpskdecmod;
        remod2 = bpskmod;
        demod2 = bpskdemod;
    elseif m2 == 4
        umod2 = qam4decmod;
        remod2 = qam4decmod;
        demod2 = qam4demod;
    elseif m2 == 16
        umod2 = qam16decmod;
        remod2 = qam16decmod;
        demod2 = qam16demod;
    elseif m2 == 64
        umod2 = qam64decmod;
        remod2 = qam64decmod;
        demod2 = qam64demod;
    end
    
    %% Choose Decoder
    if type1 == 1 || type1 == 3 || type1 == 5
        condec1 = condec12;
    elseif type1 == 7
        condec1 = condec23;
    else
        condec1 = condec34;
    end
       
    if type2 == 1 || type2 == 3 || type2 == 5
        condec2 = condec12;
    elseif type2 == 7
        condec2 = condec23;
    else
        condec2 = condec34;
    end
    
    
    
    %% Precompute some things        
    maximat = zeros(m1, m2, length(data));
    minimat = zeros(m1, m2);
    datademod1 = zeros(length(data)/log2(m1), 1);
    datademod2 = zeros(length(data)/log2(m2), 1);        
    for q = 0:(m1-1)
        for a = 0:(m2-1)
            minimat(q+1, a+1) = g1*step(umod1, q) + g2 * step(umod2, a);
        end
    end

    %% Meat of JMD
    for t = 1:length(data.')
        for q = 0:(m1-1)
            for a = 0:(m2-1)
                maximat(q+1, a+1, t) = abs(data(t) - minimat(q+1, a+1));
            end
        end
        [r, c, v] = find(maximat(:, :, t) == min(min(maximat(:,:,t))));
        datademod1(t) = r(1)-1;
        datademod2(t) = c(1)-1;
    end
    
    % Weird modulating/demodulating to get stuff to work
    dataremod1 = step(remod1, datademod1);
    dataremod2 = step(remod2, datademod2);
    datafindemod1 = step(demod1, dataremod1);
    datafindemod2 = step(demod2, dataremod2);
    
    datafindemod1 = step(convdeint, datafindemod1);
    datafindemod2 = step(convdeint, datafindemod1);
    
    datout1 = step(condec1, real(datafindemod1));
    datout2 = step(condec2, real(datafindemod2));

%%
% switch type1
%         case 1  %BPSK 1/2
%             m1 = 2;
%             switch type2
%                 case 1    
%                     m2 = 2;
%                     minimat = zeros(m1, m2, length(data));
%                     datademod1 = zeros(length(data), 1);
%                     datademod2 = zeros(length(data), 1);
%                     for t = 1:length(data.')
%                         for q = 0:(m1-1)
%                             for a = 0:(m2-1)
%                                 minimat(q+1, a+1, :) = abs(...
%                                     data(t)...
%                                     - g1*step(bpskmod, q)...
%                                     - g2*step(bpskmod, a)...
%                                     );
%                             end
%                         end
%                         [r, c, v] = find(minimat(:, :, t) ...
%                             == min(min(minimat(:,:,t))));
%                         datademod1(t) = r(1)-1;
%                         datademod2(t) = c(1)-1;
%                     end
%                     datout1 = step(condec12, real(datademod1));
%                     datout2 = step(condec12, real(datademod2));
%             end
%             
%         case 2  %BPSK 3/4
%             m1 = 2;
%             switch type2
%                 case 2
%                     m2 = 2;
%                                         minimat = zeros(m1, m2, length(data));
%                     datademod1 = zeros(length(data), 1);
%                     datademod2 = zeros(length(data), 1);
%                     for t = 1:length(data.')
%                         for q = 0:(m1-1)
%                             for a = 0:(m2-1)
%                                 minimat(q+1, a+1, :) = abs(...
%                                     data(t)...
%                                     - g1*step(bpskmod, q)...
%                                     - g2*step(bpskmod, a)...
%                                     );
%                             end
%                         end
%                         [r, c, v] = find(minimat(:, :, t) ...
%                             == min(min(minimat(:,:,t))));
%                         datademod1(t) = r(1)-1;
%                         datademod2(t) = c(1)-1;
%                     end
%                     datout1 = step(condec34, real(datademod1));
%                     datout2 = step(condec34, real(datademod2));
%         case 3  %QPSK 1/2
% 
%         case 4  %QPSK 3/4
% 
%         case 5  %16QAM 1/2
%            
%         case 6  %16QAM 3/4
%             
%         case 7  %64QAM 2/3
%             
%         case 8  %64QAM 3/4
%             
%     end
%     datout = datadec;
end