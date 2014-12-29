%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All rights reserved by Krishna Pillai, http://www.dsplog.com
% The file may not be re-distributed without explicit authorization
% from Krishna Pillai.
% Checked for proper operation with Octave Version 3.0.0
% Author        : Krishna Pillai
% Email         : krishna@dsplog.com
% Version       : 1.0
% Date          : 14th December 2008
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script for computing the BER for BPSK modulation in a
% Rayleigh fading channel with 2 Tx, 2Rx MIMO channel 
% Maximum Likelihood equalization

clear
N = 10^6; % number of bits or symbols
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

    % Maximum Likelihood Receiver
    % ----------------------------
    % if [s1 s2 ] = [+1,+1 ]
    sHat1 = [1 1];	
    sHat1 = repmat(sHat1,[1 ,N/2]);
    sHat1Mod = kron(sHat1,ones(nRx,1));	
    sHat1Mod = reshape(sHat1Mod,[nRx,nTx,N/nTx]);	
    zHat1 = squeeze(sum(h.*sHat1Mod,2)) ;
    J11 = sum(abs(y - zHat1),1);
    
    % if [s1 s2 ] = [+1,-1 ]
    sHat2 = [1 -1];	
    sHat2 = repmat(sHat2,[1 ,N/2]);
    sHat2Mod = kron(sHat2,ones(nRx,1));	
    sHat2Mod = reshape(sHat2Mod,[nRx,nTx,N/nTx]);	
    zHat2 = squeeze(sum(h.*sHat2Mod,2)) ;
    J10 = sum(abs(y - zHat2),1);
    
    % if [s1 s2 ] = [-1,+1 ]
    sHat3 = [-1 1];	
    sHat3 = repmat(sHat3,[1 ,N/2]);
    sHat3Mod = kron(sHat3,ones(nRx,1));	
    sHat3Mod = reshape(sHat3Mod,[nRx,nTx,N/nTx]);	
    zHat3 = squeeze(sum(h.*sHat3Mod,2)) ;
    J01 = sum(abs(y - zHat3),1);
    
    % if [s1 s2 ] = [-1,-1 ]
    sHat4 = [-1 -1];	
    sHat4 = repmat(sHat4,[1 ,N/2]);
    sHat4Mod = kron(sHat4,ones(nRx,1));	
    sHat4Mod = reshape(sHat4Mod,[nRx,nTx,N/nTx]);	
    zHat4 = squeeze(sum(h.*sHat4Mod,2)) ;
    J00 = sum(abs(y - zHat4),1);
    
    % finding the minimum from the four alphabet combinations 
    rVec = [J11;J10;J01;J00];
    [jj dd] = min(rVec,[],1);

    % mapping the minima to bits
    ref = [1 1; 1 0; 0 1; 0 0 ];
    ipHat = zeros(1,N);
    ipHat(1:2:end) = ref(dd,1);
    ipHat(2:2:end) = ref(dd,2);

    % counting the errors
    nErr(ii) = size(find([ip- ipHat]),2);

end

simBer = nErr/N; % simulated ber
EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbN0Lin).^(-0.5)); 
p = 1/2 - 1/2*(1+1./EbN0Lin).^(-1/2);
theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p)); 

close all
figure
semilogy(Eb_N0_dB,theoryBer_nRx1,'bp-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,theoryBerMRC_nRx2,'kd-','LineWidth',2);
semilogy(Eb_N0_dB,simBer,'mo-','LineWidth',2);
axis([0 25 10^-5 0.5])
grid on
legend('theory (nTx=1,nRx=1)', 'theory (nTx=1,nRx=2, MRC)', 'sim (nTx=2, nRx=2, ML)');
xlabel('Average Eb/No,dB');
ylabel('Bit Error Rate');
title('BER for BPSK modulation with 2x2 MIMO and ML equalizer (Rayleigh channel)');