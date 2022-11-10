% Read in Audio File
clear all;
[x,Fs] = audioread('Repeated C.aif');
t = 0:1/Fs:1-1/Fs;

% reshape to ensure single array 
[m,n] = size(x);
if n~=1
    x = reshape(x,[m*n,1]);
end

%Shortened x for debugging purposes
LX = length(x);
a = ceil(0.2*LX);
x = x(1:a);

% Add 0.25 ms Initial Silence
M = (0.25*Fs) %0.25 ms # samples
add_IS = zeros(M,1);

x = vertcat(add_IS,x);

% Plot the Noiseless and Noisy Audio Files 
Without Noise
figure();

% time domain waveform 
subplot(3,1,1, 'align')
plot(x)
title('Time Domain Waveform of Noiseless Signal x')
xlabel('Sample Number'); ylabel('Amplitude')

% FFT of audio file 
xf = fft(x);
subplot(3,1,2, 'align')
plot(abs(xf))
title('FFT of Input Noiseless Signal x')
xlabel('Sample Number'); ylabel('Amplitude')
% Spectrogram 
subplot(3,1,3, 'align')
pspectrum(x,Fs,'spectrogram','TimeResolution',0.1)
title('Spectogram of Noiseless Signal x') 
% save value into an array 
x_spec = pspectrum(x,Fs,'spectrogram','TimeResolution',0.1);

% Generate Noise Across SNR Values 
SNR = 0:0.5:20; % signal strength in dBW
for ss = 1: length(SNR)
    sigw = 1./SNR;
    N = length(x);

    y = awgn(x,SNR(ss),'measured');
    Noise = y - x;
With Noise
figure(); 
% time domain waveform 
subplot(3,1,1, 'align')
plot(y)
title('Time Domain Waveform - Noisy Signal y');
xlabel('Sample Number'); ylabel('Amplitude');

% FFT of noisy signal 
xf_n = fft(y);

subplot(3,1,2, 'align')
plot(abs(xf_n))
title('FFT of Noisy Signal y')
xlabel('Sample Number'); ylabel('Amplitude')

% spectrogram of noisy signal 
subplot(3,1,3, 'align')
pspectrum(y,Fs,'spectrogram','TimeResolution',0.1)
title('Spectogram of Noisy Signal y')
% save value into array 
y_spec = pspectrum(y,Fs,'spectrogram','TimeResolution',0.1);

% STFT
% Parameters

wind = hamming(M);
n = 1:length(x);
NN = 128;

Ts = 1/Fs;
Td = Ts*length(x)-1;
tt = linspace(0,Td,N); 
ff = linspace(0,Fs,NN);

[Sx,tx,fx] = tfrstft(x,n,NN,wind);
[Sy,ty,fy] = tfrstft(y,n,NN,wind);
[SNoise,tnoise,ynoise] = tfrstft(Noise,n,NN,wind);

figure();

% without noise 
subplot(2,1,1)
stft(x,Fs)
title('STFT of Original Signal: x')

% with noise 
subplot(2,1,2)
stft(y,Fs)
title('STFT of Noisy Signal y')

% Adaptive STFT Implementation
% notes on outcome of cqt function
%     ccfs = x_cqt
%     cvf = bandpass center frequencies
%     g = Gabor frames
%     fshifts = frequency shifts in DFT bins
%     fintervals = frequency intervals corresponding the rows of cfs
%     The kth element of fshifts is the frequency shift in DFT bins between the 
%     ((k-1) mod N) and (k mod N) element of fintervals with k = 0,1,2,...,N-1 
%     where N is the number of frequency shifts. Because MATLAB® indexes from 1,
%     fshifts(1) contains the frequency shift between fintervals{end} and 
%     fintervals{1}, fshifts(2) contains the frequency shift between fintervals{1}
%     and fintervals{2}, and so on.
%     bw = returns the bandwidth, bw, in DFT bins of the frequency intervals, fintervals
% without noise
total_BW = Fs;

L_x = length(x);

% Take fft and find norm and RMS values 
x_FFT = fft(x)/L_x;
df_x = Fs/L_x;
f_norm_x = abs(x_FFT).^2/df_x;
f_RMS_x = sqrt(df_x * sum(abs(f_norm_x)));

% store information on mean and std_dev for use later 
mean_x = sum(f_norm_x*f_RMS_x);
temp_x = (f_norm_x-mean_x).^2;
std_dev_x = sqrt((1/(L_x-1))*sum(temp_x));

% execute constant Q transform function 
[cfs_x,f_x,g_x,fshifts_x,fintervals_x,bw_x] = cqt(x);
% with noise 
total_BW = Fs;

L_y = length(y);

% Take fft and find norm and RMS values 
y_FFT = fft(y)/L_y;
df_y = Fs/L_y;
f_norm_y = abs(y_FFT).^2/df_y;
f_RMS_y = sqrt(df_y * sum(abs(f_norm_y)));

% store information on mean and std_dev for use later 
mean_y = sum(f_norm_x*f_RMS_x);
temp_y = (f_norm_x-mean_y).^2;
std_dev_y = sqrt((1/(L_y-1))*sum(temp_y));

% execute constant Q transform function 
[cfs_y,f_y,g_y,fshifts_y,fintervals_y,bw_y] = cqt(y);
% Plot results 
figure();

subplot(2,1,1)
cqt(x,'SamplingFrequency',Fs)
title('Adaptive STFT Using Constant Q Transform: Original Signal x')

subplot(2,1,2)
cqt(y,'SamplingFrequency',Fs)
title('Adaptive STFT Using Constant Q Transform: Noisy Signal y')

Spectral Subtraction on TFR 1: STFT 
Spectral Subtraction
Empty Speech is Initial 0.25 ms silence (11025 samples in x)
% STFT: Y(w) = S(w) + D(w)
Y_abs = (abs(Sy)).^2;
S_abs = (abs(Sx)).^2;

% First M columns of Sy are silence

YSP = [];
for col_index = 1:M
    YSP_temp = Sy(:,col_index);
    YSP = [YSP YSP_temp];
end

Dhat_abs = mean(abs(YSP)).^2;
[a b] = size(Y_abs);
repetitions = ceil(b/M);
Dhat_abs = repmat(Dhat_abs,a,repetitions);
[g h] = size(Dhat_abs);
rid_columns = h-b
rid_elements = rid_columns*g
Dhat_abs = Dhat_abs(1:((g*h)-rid_elements));
Dhat_abs = reshape(Dhat_abs,a,b);
% Real D
D_abs = Y_abs - S_abs;
D_abs2 = (abs(SNoise)).^2;

H2w = 1-(Dhat_abs)./Y_abs;
Shat_abs = H2w.*Y_abs;
b_temp = max(0,H2w);
Hw = sqrt(b_temp);
 
Sxhat = Hw.*Sy; 

% Recover Signal and Plot in Time Domain
xhat_STFT = tfristft(Sxhat,tx,wind);

figure();
subplot(3,1,1) 
plot(x)
title('Original Noiseless Signal x')
xlabel('time'); ylabel('amplitude');
subplot(3,1,2) 
plot(xhat_STFT)
title('Estimated Signal x from Spectral Subtraction - STFT')
xlabel('time'); ylabel('amplitude');

subplot(3,1,3)
plot(y)
title('Time Domain Waveform - Noisy Signal y');
xlabel('Sample Number'); ylabel('Amplitude');
Plot in Frequency Domain
% without noise 
figure();
subplot(2,1,1)
stft(y,Fs)
title('STFT of Noisy Signal y ')
% with noise 
subplot(2,1,2)
stft(xhat_STFT,Fs)
% title('Estimated STFT of Original Signal')
title('Spectral Subtraction Estimation of Noiseless Signal x Using STFT Data')

Spectral Subtraction on TFR 2: Adaptive STFT 
Empty Speech is Initial 0.25 ms silence (11025 samples in x)
% STFT: Y(w) = S(w) + D(w)
Y_abs_ASTFT = (abs(cfs_y)).^2;
S_abs_ASTFT = (abs(cfs_x)).^2;

% First M columns of Sy are silence

YSP2 = [];
for col_index2 = 1:M
    YSP_temp2 = Sy(:,col_index2);
    YSP2 = [YSP2 YSP_temp2];
end

Dhat2_abs = mean(abs(YSP2)).^2;
[a2, b2] = size(Y_abs_ASTFT);
repetitions2 = ceil(b2/M);
Dhat2_abs = repmat(Dhat2_abs,a2,repetitions2);
[g2, h2] = size(Dhat2_abs);
rid_columns2 = h2-b2
rid_elements2 = rid_columns2*g2
Dhat2_abs = Dhat2_abs(1:((g2*h2)-rid_elements2));
Dhat2_abs = reshape(Dhat2_abs,a2,b2);
% Real D
D_abs_ASTFT = Y_abs_ASTFT - S_abs_ASTFT;
D_abs2_ASTFT = (abs(SNoise)).^2;

H2w_ASTFT = 1-(Dhat2_abs)./Y_abs_ASTFT;
Shat_abs_ASTFT = H2w_ASTFT.*Y_abs_ASTFT;
b_temp2 = max(0,H2w_ASTFT);
Hw2 = sqrt(b_temp2);
 
Sxhat2 = Hw2.*cfs_y; 

Recover Signal and Plot in Time Domain
wind2=220687;
xhat_ASTFT = icqt(Sxhat2,g_x,fshifts_x);

figure();
subplot(3,1,1) 
plot(x)
title('Original Noiseless Signal x')
xlabel('time'); ylabel('amplitude');

subplot(3,1,2) 
plot(xhat_ASTFT)
title('Estimated Signal x From Spectral Subtraction - Adaptive STFT')
xlabel('time'); ylabel('amplitude');

subplot(3,1,3) 
plot(y)
title('Time Domain Waveform - Noisy Signal y');
xlabel('Sample Number'); ylabel('Amplitude');
Plot in Frequency Domain
figure(); 
% without noise 
subplot(2,1,1)
stft(y,Fs)
title('STFT of Noisy Signal y')

% with noise 
subplot(2,1,2)
stft(xhat_ASTFT,Fs)
%title('Estimated STFT of Original Signal Using Spectral Subtraction on Adaptive STFT')
title('Spectral Subtraction Estimation of Noiseless Signal x Using Adaptive STFT Data')
Efficiency Analysis - STFT  
    errors = 0;
    L = length(xhat_STFT);

    % x_actual = x(1:L); % disregard size decrease from STFT
    errors=0;
    range = 0.05;
    diff_STFT = abs(x - xhat_STFT);
    %range = (max(diff_STFT)-min(diff_STFT))/length(diff_STFT);
    for i = 1 : L
        if diff_STFT(i) > range
            errors = errors+1;
        end
    end

    percent(ss) = (errors/L)*100;

end


semilogy(SNR,percent, 'linewidth',2);
grid on
xlabel('SNR (dB)');
ylabel('BER');

title('BER Vs SNR graph For STFT')
Efficiency Analysis - Adaptive STFT  
    errors = 0;
    L = length(xhat_ASTFT);

    % x_actual = x(1:L); % disregard size decrease from STFT
    errors=0;
    range = 0.05;
    diff_STFT = abs(x - xhat_ASTFT);
    %range = (max(diff_STFT)-min(diff_STFT))/length(diff_STFT);
    for i = 1 : L
        if diff_STFT(i) > range
            errors = errors+1;
        end
    end

    percentA(ss) = (errors/L)*100;

end

semilogy(SNR,percentA, 'linewidth',2);
grid on
hold on
semilogy(SNR,percent(1:39), 'linewidth',2);
xlabel('SNR (dB)');
ylabel('BER');
xlim([1 20])
ylim([0 179])
title('BER Vs SNR graph For STFT') 
