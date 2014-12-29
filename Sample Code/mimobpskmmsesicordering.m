%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All rights reserved by Krishna Pillai, http://www.dsplog.com
% The file may not be re-distributed without explicit authorization
% from Krishna Pillai.
% Checked for proper operation with Octave Version 3.0.0
% Author        : Krishna Pillai
% Email         : krishna@dsplog.com
% Version       : 1.0
% Date          : 06 December 2008
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script for computing the BER for BPSK modulation in a
% Rayleigh fading channel with 2 Tx, 2Rx MIMO channel 
% Minimum Mean Square Error Equalization with Successive Interference 
% Cancellation (ZF-SIC) with optimal ordering

clear
N =  10^6; % number of bits or symbols
Eb_N0_dB = [0:25]; % multiple Eb/N0 values
nTx = 2;
nRx = 2;
for ii = 1:length(Eb_N0_dB)

    % Transmitter
    ip = rand(1,N)>0.5; % generating 0,1 with equal probability
    s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 0
    sMod = kron(s,ones(nRx,1)); % 
    sMod = reshape(sMod,[nRx,nTx,N/nTx]); % grouping in [nRx,nTx,N/NTx ] matrix
    
    h = 1/sqrt(2)*[randn(nRx,nTx,N/nTx) + j*randn(nRx,nTx,N/nTx)]; % Rayleigh channel
    n = 1/sqrt(2)*[randn(nRx,N/nTx) + j*randn(nRx,N/nTx)]; % white gaussian noise, 0dB variance

    % Channel and noise Noise addition
    y = squeeze(sum(h.*sMod,2)) + 10^(-Eb_N0_dB(ii)/20)*n;

    % Receiver
    % ----------
    % Forming the MMSE equalization matrix W = inv(H^H*H + sigma^2*I)*H^H
    % H^H*H is of dimension [nTx x nTx]. In this case [2 x 2] 
    % Inverse of a [2x2] matrix [a b; c d] = 1/(ad-bc)[d -b;-c a]
    hCof = zeros(2,2,N/nTx)  ; 
    hCof(1,1,:) =  sum(h(:,2,:).*conj(h(:,2,:)),1) + 0*10^(-Eb_N0_dB(ii)/10);  % d term
    hCof(2,2,:) =  sum(h(:,1,:).*conj(h(:,1,:)),1) + 0*10^(-Eb_N0_dB(ii)/10);  % a term
    hCof(2,1,:) = -sum(h(:,2,:).*conj(h(:,1,:)),1); % c term
    hCof(1,2,:) = -sum(h(:,1,:).*conj(h(:,2,:)),1); % b term
    
    for kk = 1:2
    
       if kk == 1
          sortIdx = [];
    	  hCof(1,1,:) =  sum(h(:,2,:).*conj(h(:,2,:)),1) + 10^(-Eb_N0_dB(ii)/10);  % d term
          hCof(2,2,:) =  sum(h(:,1,:).*conj(h(:,1,:)),1) + 10^(-Eb_N0_dB(ii)/10);  % a term
          hCof(2,1,:) = -sum(h(:,2,:).*conj(h(:,1,:)),1); % c term
          hCof(1,2,:) = -sum(h(:,1,:).*conj(h(:,2,:)),1); % b term
       elseif kk == 2
          % Sorting the equalization matrix based on the channel power on each dimension
          % since the second spatial dimension is equalized first, the channel
          % with higher power assigned to second dimension
          normSS1 = squeeze(hCof(2,2,:));
    	  normSS2 = squeeze(hCof(1,1,:));
    	  sortIdx = find(normSS2 < normSS1);
	end

   
        % sorting the H^H*H  + sigma^2*I matrix 
        hCofSort = hCof;
        if ~isempty(sortIdx)
            hCofSort(2,2,sortIdx) = hCof(1,1,sortIdx) + 10^(-Eb_N0_dB(ii)/10);;
      	    hCofSort(1,1,sortIdx) = hCof(2,2,sortIdx) + 10^(-Eb_N0_dB(ii)/10);;
    	    hCofSort(1,2,sortIdx) = hCof(2,1,sortIdx);
    	    hCofSort(2,1,sortIdx) = hCof(1,2,sortIdx);
        end
        hDen = ((hCofSort(1,1,:).*hCofSort(2,2,:)) - (hCofSort(1,2,:).*hCofSort(2,1,:))); % ad-bc term
        hDen = reshape(kron(reshape(hDen,1,N/nTx),ones(2,2)),2,2,N/nTx);  % formatting for division
        hInvSort = hCofSort./hDen; % inv(H^H*H)

        % sorting the H matrix
        hSort = h;
        if ~isempty(sortIdx)
    	    hSort(:,2,sortIdx) = h(:,1,sortIdx);
    	    hSort(:,1,sortIdx) = h(:,2,sortIdx);
        end

        % Equalization - Zero forcing
        hModSort =  reshape(conj(hSort),nRx,N); % H^H operation
    
        yModSort = kron(y,ones(1,2)); % formatting the received symbol for equalization
        yModSort = sum(hModSort.*yModSort,1); % H^H * y 
        yModSort =  kron(reshape(yModSort,2,N/nTx),ones(1,2)); % formatting
        yHatSort = sum(reshape(hInvSort,2,N).*yModSort,1); % inv(H^H*H)*H^H*y

        % receiver - hard decision decoding on second spatial dimension
        ipHat2SS = real(yHatSort(2:2:end))>0;
        ipHatMod2SS = 2*ipHat2SS-1;
        ipHatMod2SS = kron(ipHatMod2SS,ones(nRx,1));
        ipHatMod2SS = reshape(ipHatMod2SS,[nRx,1,N/nTx]);

        % new received symbol - removing the effect from second spatial dimension
        h2SS = hSort(:,2,:); % channel in the second spatial dimension
        r = y - squeeze(h2SS.*ipHatMod2SS);

        % maximal ratio combining - for symbol in the first spatial dimension
        h1SS = squeeze(hSort(:,1,:));
        yHat1SS = sum(conj(h1SS).*r,1)./sum(h1SS.*conj(h1SS),1);
        yHatSort(1:2:end) = yHat1SS;
  
        yHatSort = reshape(yHatSort,2,N/2) ;
        if ~isempty(sortIdx)
            yHatSort(:,sortIdx) = flipud(yHatSort(:,sortIdx));
        end
        yHat = reshape(yHatSort,1,N);

        % receiver - hard decision decoding
        ipHat = real(yHat)>0;

        % counting the errors
        nErr(kk,ii) = size(find([ip- ipHat]),2);

    end

end

simBer = nErr/N; % simulated ber
EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbN0Lin).^(-0.5)); 
p = 1/2 - 1/2*(1+1./EbN0Lin).^(-1/2);
theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p)); 

close all
semilogy(Eb_N0_dB,theoryBer_nRx1,'bp-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,theoryBerMRC_nRx2,'kd-','LineWidth',2);
semilogy(Eb_N0_dB,simBer(1,:),'mo-','LineWidth',2);
semilogy(Eb_N0_dB,simBer(2,:),'gp-','LineWidth',2);
axis([0 25 10^-5 0.5])
grid on
legend('theory (nTx=2,nRx=2, ZF)', 'theory (nTx=1,nRx=2, MRC)', 'sim (nTx=2, nRx=2, MMSE-SIC)','sim (nTx=2, nRx=2, MMSE-SIC-Sort)');
xlabel('Average Eb/No,dB');
ylabel('Bit Error Rate');
title('BER for BPSK modulation with 2x2 MIMO and MMSE-SIC equalizer (Rayleigh channel)');


