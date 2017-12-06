clear all;
close all;

%% CONSTANTS

% sampling information
Fs = 16000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    


%% LOAD AND STORE THE DATA FILES

S1 = load('data.mat', 'Clean');
S2 = load('data.mat', 'Data');

sig_clean = S1.Clean; % clean signal (1 channel)
sig_noisy_M = S2.Data; % raw signal (16 channels)

M = size(sig_noisy_M, 2); % number of channels

%%  COMPUTE THE DISRETE-FOURIER-TRANSFORM

S = dft(sig_clean,Fs);      % STFT of clean sig


Y = zeros([size(S) M]); %preallocate space for the 16 channels STFT

for i = 1 : M       % store the 16 STFTs in a vector of matrices
     Y(:, :, i) = dft(sig_noisy_M(:,  i),Fs);
end

%% A FIRST ESTIMATOR: SAMPLE MEAN

Y_avg = mean(Y,3); %% Average over the 16 channels (third dimension)


%% COMPUTE THE INVERSE-DISCRETE-FOURIER-TRANSFORM

y_avg = idft(Y_avg, Fs);
s = idft(S, Fs);

%% PLOT OR PlAY

% %time plot
% t = (1:length(s))*T;
% figure;
% plot(t,s, 'b');
% title('Clean signal')

% frequency plot
% Yavg = mean(S,2);
% Py = abs(Yavg);
% x = T*(1:length(Py));
% figure
% plot(x,Py)

% sound( s (10e3:100e3), Fs)

%% ERRORS

sig_clean_cut = sig_clean(1 : length(y_avg) ); % adjust the length for plotting
err = (sig_clean_cut - y_avg); % error using a SAMPLE MEAN ESTIMATOR
t = (1:(length(err)))*T;
figure
subplot(2,1,1)
plot(t, err)
title('Error clean sig - sample mean estimation');
axis([-inf inf -5 5]);


err2 = (sig_clean - sig_noisy_M(:,4)); % error due to the STFT
t2 = (1:(length(err2)))*T;
subplot(2,1,2)
plot(t2, err2)
title('Error clean sig - unprocessed noisy sig');
axis([-inf inf -5 5]);