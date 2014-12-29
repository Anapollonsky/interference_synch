%% Andrew Apollonsky
clc;
close all;
clear all;

len = 1000;

% Generate Sent Signal
sig = round(rand(1, len))*2 - 1; % Generate -1/+1 value
sigbin = sig/2 + .5; % Convert to 0/1 binary

% Generate Noise
noipwrdb = -30:.3:0;  
noipwr = 10.^(fliplr(noipwrdb)/10);
noi = repmat(sqrt(noipwr.'), [1 len]) .* randn(length(noipwrdb), len);

% Generate Interference  
int = round(rand(1, len))*2 - 1;  
sirdb = [-5, 0, 5];
sir = 10.^(sirdb/10);

% Generate other values
snrdb = -fliplr(noipwrdb);


% Properly scale and such 
recsig = repmat(sig, [length(noipwrdb) 1 length(sir)]);
recint = shiftdim(repmat(int.' * sir, [1 1 length(noipwrdb)]), 2);
recnoi = repmat(noi, [1 1 length(sir)]);

rec = recsig + recint + recnoi;

%% Interference-Blind
decir = interp1([-1000 1000],[-1 1],rec,'nearest'); %Find nearest -1/+1 value
decib = decir/2 + .5; % Convert to 0/1 binary


syeray = zeros(length(sirdb), length(noipwrdb));

for k = 1:length(sirdb) % For every SIR
    for m = 1:length(noipwrdb)  % For every SNR
        [unneeded, syeray(k,m)] = biterr(sigbin, decib(m,:,k));
    end
end
%% Plot
figure;
hold all;
for k = 1:length(sirdb) % For every SIR
    plot(snrdb, syeray(k,:));
end
xlabel('SNR');
ylabel('SER');
legend('5', '0', '-5');
set(gca,'yscale','log');
hold off;
    
% Rows: noise pwr. 
% Columns: Time
% Depth: SIR