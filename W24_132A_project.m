M=4;         % QPSK
m=2;         % represent each symbol with 2 bits
Es=1;
Nbits = 1e6; % number of bits
Rs = 1; % symbol rate
Fs = 4; % sampling rate
sps = Fs/Rs; % samples per symbol
EsN0dB = -3:1:20;
% EsN0dB=10;
EbN0dB = EsN0dB-3; % Es=kEb, k=2 for QPSK, 10log10(2)~= 3
N0 = 10.^(-EsN0dB/10);
% b = rcosdesign(beta,span,sps,shape)
% hrrc = rcosdesign(0.2,12,sps,"normal");
hrrc = rrc((-6:1/sps:6)',0.2,1);

% Tx
% 1a: Generate a sufficient long vector of bits
bitvec = randi([0,1],[Nbits,1]);

% 1b: Generate a Gray-coded QPSK vector of complex symbols
% 1/sqrt(2) normalize symbol energy to make average symbol energy Es=1
symvec = complex(2*bitvec(1:2:end)-1,2*bitvec(2:2:end)-1)/sqrt(2);

%% AWGN channel
% 1c: Generate the "complex baseband sampled waveform" by applying SRRC pulse shaping. Upsample by at least a factor 2. 
txvec = conv(upsample(symvec,sps),hrrc,'same');

% 1d(i): AWGN channel, Add noise to the modulated waveform to match a desire
% SNR from -3dB to 20dB
Gray_mapping=[1 1; -1 1; -1 -1; 1 -1];
QPSK_modulation=[1 1; -1 1; -1 -1; 1 -1]*(1/sqrt(2));
BER=zeros(1,length(EsN0dB));
for i=1:length(EsN0dB)
    % only noise will alter with SNR
    sig2 = N0(i)*Fs/2;
    % For each SNR, we have 2*1e6 trasmitted points
    noisevec = sqrt(sig2)*complex(randn(size(txvec)),randn(size(txvec)));
    % Rx
    rxvec = txvec + noisevec; % length=(2e6x1)
    
    % 1d(ii): Pass the received signal rxvec through a match filter, hrrc(T-t),
    % and sample at one sample per symbol

    % y(t)=rxvex(i,:)*hrrc
    matched_filter=conv(rxvec,hrrc,'same');
    % downsample y(t)
    matched_filter_sample=matched_filter(1:sps:end);

    % 1d(iii): Demodulate
    decoded_bits=zeros(1,length(bitvec));
    for j=1:length(matched_filter_sample)
        received_point=[real(matched_filter_sample(j)) imag(matched_filter_sample(j))];
        %{ 
        % method i
        % Loop through each point in the constellation to find the closest one
        distance=zeros(1,M);
        for k=1:M
            % norm(b-a): || b-a ||
            % distance(1,k)=sqrt((QPSK_modulation(k,1)-received_point(1))^2+(QPSK_modulation(k,2)-received_point(2))^2);
            distance(1,k)=norm(Gray_mapping(k,:)-received_point);
        end
        [min_distance, index]=min(distance);
        decoded_bits(1,(2*j-1):2*j)=(Gray_mapping(index,:)+1)/2;
        %}
        
        %{
        % method ii
        if (received_point(1)>0) && (received_point(2)>0)
            decoded_bits(1,(2*j-1):2*j)=(Gray_mapping(1,:)+1)/2;
        elseif (received_point(1)<0) && (received_point(2)>0)
            decoded_bits(1,(2*j-1):2*j)=(Gray_mapping(2,:)+1)/2;
        elseif (received_point(1)<0) && (received_point(2)<0)
            decoded_bits(1,(2*j-1):2*j)=(Gray_mapping(3,:)+1)/2;
        elseif (received_point(1)>0) && (received_point(2)<0)
            decoded_bits(1,(2*j-1):2*j)=(Gray_mapping(4,:)+1)/2;
        end
        %}
        
    end
    [syms,decoded_bits]=qpskdetect(matched_filter_sample);
    BER(i)=sum(xor(bitvec,decoded_bits))/Nbits;
end

% theoretical P(e) for gray-mapping QPSK
EsN0dB2 = -3:0.1:20;
EbN0dB2 = EsN0dB2 - 3;
N02 = 10.^(-EsN0dB2/10);
BER_th=qfunc(sqrt(Es./N02)).*(1-0.5.*qfunc(sqrt(Es./N02)));

clf;    % Clears the current figure window
figure; % Opens a new figure window
scatter(EbN0dB,BER,'filled')
set(gca,'yscale','log')
hold on
plot(EbN0dB2,BER_th,'LineWidth',1.3)
xlim([-6 20])
ylim([10e-7 1])
grid on
hold off
xlabel('E_b/N_0 (dB)', 'Interpreter', 'tex');
ylabel('Bit Error Rate (BER)');
title('BER vs. Eb/N0(dB) for QPSK through AWGN channel')
legend('Emperical BER','Theoretical BER')




% Plot BER v.s. Eb/N0(dB)
% documentation: readme
% testing from different aspects, ex: QPSK modulation library (comm.QPSKModulator)