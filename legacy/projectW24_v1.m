Nbits = 1e6; % number of bits
Rs = 1; % symbol rate
Fs = 4; % sampling rate
sps = Fs/Rs; % samples per symbol
EsN0dB = 10;
EbN0dB = EsN0dB-3;
N0 = 10.^(-EsN0dB/10);
hrrc = rrc((-6:1/sps:6)',0.2,1);
% Tx
bitvec = randi([0,1],[Nbits,1]);
symvec = complex(2*bitvec(1:2:end)-1,2*bitvec(2:2:end)-1)/sqrt(2);
%% AWGN channel
% sample-level simulation
txvec = conv(upsample(symvec,sps),hrrc,'same');
sig2 = N0*Fs/2;
noisevec = sqrt(sig2)*complex(randn(size(txvec)),randn(size(txvec)));
% Rx
rxvec = txvec + noisevec;
% MF
mfout = ...; % matched filter rxvec
mfsamp = ...; % sample mfout
[detsymvec,detbitvec] = qpskdetect(mfsamp);
ser1 = sum(detsymvec~=symvec)/(Nbits/2);
ber1 = sum(xor(detbitvec,bitvec))/Nbits;
% symbol level
symnoise = sqrt(N0/2)*complex(randn(size(symvec)),randn(size(symvec)));
rxsymvec = symvec + symnoise;
[detsymvec2,detbitvec2] = qpskdetect(rxsymvec);
ser2 = sum(detsymvec2~=symvec)/(Nbits/2);
ber2 = sum(xor(detbitvec2,bitvec))/Nbits;
% compare sample-level to symbol-level sims
% repeat for different SNRs
%% Flat fading
% sample-level simulation
txvec = conv(upsample(symvec,sps),hrrc,'same');
sig2 = N0*Fs/2;
sigfad = 1;
tmps = 0;
tmpb = 0;
nruns = 100;
for k = 1:nruns
    noisevec = sqrt(sig2)*complex(randn(size(txvec)),randn(size(txvec)));
    fadingalpha = sigfad*complex(randn(1),rand(1));
    rxvecfad = fadingalpha * txvec + noisevec;
    mfout = conv(rxvecfad*exp(-1i*angle(fadingalpha)),hrrc,'same')/Fs; % phase compensation
    mfsamp = mfout(1:sps:end);
    [detsymvec,detbitvec] = qpskdetect(mfsamp);
    tmps = tmps + sum(detsymvec~=symvec)/(Nbits/2);
    tmpb = tmpb + sum(xor(detbitvec,bitvec))/Nbits;
end
ser3 = tmps/nruns;
ber3 = tmpb/nruns;
% repeat to average fading
% repeat for different SNRs
% symbol-level
% do this
%% 2x2 MIMO
Nt = 2; Nr = 2; % number of tx and rx antennas
sigfad = 1;
H = complex(randn(Nr,Nt),randn(Nr,Nt)); % channel matrix
txsyms = reshape(symvec,floor([Nt,length(symvec)/Nt]));
% perfect CSI at Tx and Rx - MIMO for throughput
[u,s,v] = svd(H);
xtilde = ... ; % precode txsyms;
noisemat = sqrt(N0)*complex(randn(size(txsyms)),randn(size(txsyms)));
ytilde = H*xtilde + noisemat;
y = ...; % postprocess ytilde
rxvec = y(:);
[detsymvec,detbitvec] = qpskdetect(rxvec);
ser4 = sum(detsymvec~=symvec)/(Nbits/2);
ber4 = sum(xor(detbitvec,bitvec))/Nbits;
% perfect CSI at Tx and Rx - MIMO for diversity
txsyms = repmat(symvec,[1,2]).'; % note the difference with line 51!
xtilde = ... ; % precode txsyms;
noisemat = sqrt(N0)*complex(randn(size(txsyms)),randn(size(txsyms)));
ytilde = H*xtilde + noisemat;
y = ...; % postprocess ytilde
rxvec = sum(y).';
[detsymvec,detbitvec] = qpskdetect(rxvec);
ser5 = sum(detsymvec~=symvec)/(Nbits/2);
ber5 = sum(xor(detbitvec,bitvec))/Nbits;
%% perfect CSI at Rx only
% ZF - throughput
txsyms = reshape(symvec,floor([Nt,length(symvec)/Nt]));
noisemat = sqrt(N0)*complex(randn(size(txsyms)),randn(size(txsyms)));
y = H*txsyms + noisemat;
xhatzf = ...; % apply ZF channel estimate
rxvec = xhatzf(:);
[detsymvec,detbitvec] = qpskdetect(rxvec);
ser6 = sum(detsymvec~=symvec)/(Nbits/2);
ber6 = sum(xor(detbitvec,bitvec))/Nbits;
% LMMSE - throughput
xhatlmmse = ...; % apply LMMSE channel estimate
rxvec = xhatlmmse(:);
[detsymvec,detbitvec] = qpskdetect(rxvec);
ser7 = sum(detsymvec~=symvec)/(Nbits/2);
ber7 = sum(xor(detbitvec,bitvec))/Nbits;
% ZF - diversity
% do this
% LMMSE - diversity
% do this
%% Estimate channel
% use pilots (QPSK symbols, but different from info symbols)
% detect using estimated channel
% do this
%% Alamouti (Nt = 2)
txsyms = zeros(2,length(symvec));
for ii = 1:length(symvec)/2
    txsyms(:,(ii-1)*2+1) = symvec((ii-1)*2+(1:2));
    txsyms(:,(ii-1)*2+2) = [-conj(symvec((ii-1)*2+2));conj(symvec((ii-1)*2+1))];
end
noisemat = sqrt(N0)*complex(randn(size(txsyms)),randn(size(txsyms)));
y = H*txsyms + noisemat;
yala = zeros(4,length(symvec)/2);
for ii = 1:length(symvec)/2
    tmp = y(:,(ii-1)*2+(1:2));
    yala(:,ii) = [tmp(1,1);conj(tmp(1,2));tmp(2,1);conj(tmp(2,2))];
end
Hala = [H(1,1),H(1,2);conj(H(1,2)),-conj(H(1,1));H(2,1),H(2,2);conj(H(2,2)),-conj(H(2,1))];
tmp1 = Hala'*yala;
rxvec = reshape(tmp1,[1,length(symvec)]).';
[detsymvec,detbitvec] = qpskdetect(rxvec);
ser8 = sum(detsymvec~=symvec)/(Nbits/2);
ber8 = sum(xor(detbitvec,bitvec))/Nbits;
% try Nr = 1, 2, 3
