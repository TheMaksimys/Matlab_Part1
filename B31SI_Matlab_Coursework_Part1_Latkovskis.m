%% This is the code for the MATLAB assignment Part 1 for B31SI course
% The program emulates the wireless transmission of the QPSK-coded
% bitstream over the Rayleigh channel with AWGN. The emulation result is
% the Bit-Error-Rate versus the SNR of the recieved signal plot.

clear all
clf

NumberOfBits=500000; % Length of the transmitted bit sequence over the channel
BitRate = 541600; % Bits/second
fs = 270800; % Sampling frequency
Ts = NumberOfBits/BitRate-1/fs; % Message transmission period

%% Bit source
Vhigh=1; % Voltage level for HIGH bit state (Here normalized to 1 Volt)
Vlow=0;  % Voltage level for LOW bit state

BinSeq=randi([Vlow Vhigh],1,NumberOfBits); % Generation of random binary bit sequence
[BinSeqmin, BinSeqmax]=bounds(BinSeq); % Function to check upper and lower boundaries of the generated sequence 

%% QPSK Modulation (Grey coding)
BinOfs=(Vhigh-BinSeq*2); % Converting binary states  to [-a;a] format
BinSeqIQ=zeros(2,NumberOfBits/2); % Creating an empty matrix for I and Q component allocation

BinSeqIQ(1,:)=BinOfs(:,1:2:end);  % Deternmining the In-phase component 
BinSeqIQ(2,:)=BinOfs(:,2:2:end);  % Deternmining the Quadrature component

S=zeros(1,NumberOfBits/2); % Signal to be transmitted through the channel
Theta=3.14/4; % Angle on the Constellation diagram
S(1,:)=1j*sin(Theta).*BinSeqIQ(1,:)+cos(Theta).*BinSeqIQ(2,:); % Gray coding of the input bitstream

%% Rayleigh Fading Channel

b=0.5; % Channel variance
fm=91; % Maximum Doppler frequency
N1=9; % Number of sinusoids for In-phase component
N2=10; % Number of sinusoids for Quadrature component

t=[0:1/fs:Ts]; % Time divisions for the channel simulation

c1=zeros([1 N1]); % Generating zero vector for the In-phase sinusoinds amplitude
c2=zeros([1 N2]); % Generating zero vector for the Quadrature sinusoinds amplitude
f1=zeros([1 N1]); % Generating zero vector for the In-phase sinusoinds frequency offset
f2=zeros([1 N2]); % Generating zero vector for the Quadrature sinusoinds frequency offset
theta1=zeros([1 N1]); % Generating zero vector for the In-phase sinusoid phase offset
theta2=zeros([1 N2]); % Generating zero vector for the Quadrature sinusoid phase offset

% Computing parameters for the In-phase sinusoids
for n = 1:1:N1
    c1(1,n) = sqrt(2*b/N1); % Amplitude
    theta1(1,n)=2*pi*(n/(N1+1)); % Phase
    f1(1,n)=fm*sin(((pi/2)/N1)*(n-0.5)); % Frequency
end

% Computing parameters for the Quadrature sinusoids
for n = 1:1:N2
    c2(1,n) = sqrt(2*b/N2); % Amplitude
    theta2(1,n)=2*pi*(n/(N2+1)); % Phase
    f2(1,n)=fm*sin(((pi/2)/N2)*(n-0.5)); % Frequency
end

g1_n=zeros(N1, numel(t)); % Generating zero matrix for In-phase sinusoids
g2_n=zeros(N2, numel(t)); % Generating zero matrix for Quadrature sinusoids
 
% Generating In-phase and Quadrature sinusoids
for i=1:1:numel(t)
    g1_n(:,i)=c1.*cos(2*pi.*f1*t(1,i)+theta1);
    g2_n(:,i)=c2.*cos(2*pi.*f2*t(1,i)+theta2);
end

g1=sum(g1_n,1); % Sum of In-phase sinusoids
g2=sum(g2_n,1); % Sum of Quadrature sinusoids

G(1,:)=g1(1,:)+1j*g2(1,:); % Computing Rayleigh flat fading channel

%% AWGN

SNR = [0:2:20]; % SNR of the signal in dB
SNR_lin=10.^(SNR/10); % SNR in linear scale
k=2; % Number of bits in one symbol

BERsimulated=zeros(1,numel(SNR)); % Creating zero vector to store simulated BER

for m=1:1:numel(SNR) % Performing simulations for different SNRs

    NormRealAWGN=randn([1 numel(S)]); % Generating real component for AWGN simulation
    NormImAWGN=randn([1 numel(S)]); % Generating imaginary component for AWGN simulation

    NormAWGN=(NormRealAWGN+1j*NormImAWGN); % Computing normalized AWGN with 1/2 variance
    NoiseAmp=Vhigh/(2*SNR_lin(1,m)*k); % Computing the Average Noise power for the AWGN
    AWGN=sqrt(NoiseAmp)*NormAWGN; % Computing AWGN for the given SNR

    %% Emulation of the recieved signal through AWGN Rayleigh communication channel

    R=zeros(1,numel(S)); % Creating the Recieved signal matrix
    R=S.*G+AWGN; % Computing the recieved signal matrix in accordance with the system diagram

    %% Signal demodulation

    Rhat=zeros(1,numel(R)); % Creating zeros vector to store the reconstructed signal
    Rhat=conj(G).*R; % Reconstructing the transmitted signal
    Gnorm=conj(G).*G; % Computing modulus of G
    Rhat=Rhat./Gnorm; % Dividing the reconstructed signal by the modulus of G to obtain correct constellation diagram

%     % This plot is purely for debugging purposes and shall not be used during the normal run
%     figure (3) % Plotting everything on the figure 2
%     plot(real(Rhat),imag(Rhat),'m.'); % Plotting the reconstructed signal on the constellation diagram
%     grid on
%     title('Reconstructed signal constellation');

    % Listing all allowed QPSK constellation points on the IQ diagram for the Gray coding
    Constellation00=Vhigh*exp(1j*pi/4); % Represents "00" codeword
    Constellation10=Vhigh*exp(1j*3*pi/4); % Represents "01" codeword
    Constellation11=Vhigh*exp(1j*5*pi/4); % Represents "11" codeword
    Constellation01=Vhigh*exp(1j*7*pi/4); % Represents "10" codeword
 
    MLerror=zeros(4,numel(Rhat)); % Creating Maximum likelyhood matrix
    % Computing the absolute errors with respect to the allowed constellation points
    MLerror=[abs(Rhat-Constellation00);abs(Rhat-Constellation10);abs(Rhat-Constellation11);abs(Rhat-Constellation01)];

    RecSignal=zeros(2,numel(Rhat)); % Creating the recovered signal matrix
    [RecSignal(1,:), RecSignal(2,:)] = min(MLerror); % Recovering signal using ML technique
 
    BitStreamOut=zeros(1,NumberOfBits); % Creating output bitstream vector
    for i=1:1:numel(R)
        if(RecSignal(2,i)==1) % If 00 codeword is decoded, assign:
             BitStreamOut(1,2*i-1)=0; % Bit 0 to 0
             BitStreamOut(1,2*i)=0; % Bit 1 to 0
        elseif(RecSignal(2,i)==2) % If 01 codeword is decoded, assign:
             BitStreamOut(1,2*i-1)=0; % Bit 0 to 0
             BitStreamOut(1,2*i)=1; % Bit 1 to 1
         elseif(RecSignal(2,i)==3) % If 11 codeword is decoded, assign:
             BitStreamOut(1,2*i-1)=1; % Bit 0 to 1
             BitStreamOut(1,2*i)=1; % Bit 1 to 1
         elseif(RecSignal(2,i)==4) % If 10 codeword is decoded, assign:
             BitStreamOut(1,2*i-1)=1; % Bit 0 to 1 
             BitStreamOut(1,2*i)=0; % Bit 0 to 0
        end
    end
 
    %% BER computaion by comparing the bits at the bit source and the bit sink
    BERsimulated(1,m) = (sum(BinSeq~=BitStreamOut))/NumberOfBits; % Calculating the Bit Error Ratio
end

%% Couputing theoretical BER

BERtheoretical=zeros(1,numel(SNR)); % Creating zero vector to store theoretical BER
gamma=zeros(1,numel(SNR)); % Crearing zero vector to store the gamma coefficients

for m=1:1:numel(SNR)
 gamma(1,m)=2*b*SNR_lin(1,m); % Computing the gamma coefficient for the BER calculation
 BERtheoretical(1,m)=0.5*(1-sqrt(gamma(1,m)/(1+gamma(1,m)))); % Calculating the theoretical BER for Gray coded QPSK signal in the slow fading Rayleigh channel
end

%% Plotting simulated and analytical results

% Plotting Analytical and Simulated BER in the semilog scale
figure (1) % Pick figure 1
semilogy(SNR,BERsimulated); % Plotting simulated BER with the logaritmic Y axis
hold on
semilogy(SNR,BERtheoretical); % Plotting theoretical BER with the logaritmic Y axis
hold off
grid on;
legend('simulation','theory'); % Labelling plots
xlabel('Eb/No (dB)'); % Labelling X axis
ylabel('Bit Error Rate'); % Labelling Y axis
title('BER for Gray QPSK over Rayleigh channel'); % Labelling figure

%% End of the document