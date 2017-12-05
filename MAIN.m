clear all;

%% CONSTANTS

% sampling information
Fs = 16000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    


%% LOAD AND STORE THE DATA FILES

S1 = load('data.mat', 'Clean');
S2 = load('data.mat', 'Data');

sig_clean = S1.Clean; % clean signal (1 channel)
sig_noisy = S2.Data(:,1); % raw signal (16 channels)

%%  COMPUTE THE DISRETE-FOURIER-TRANSFORM

Y = dft(sig_noisy,Fs);
S = dft(sig_clean,Fs);


%% COMPUTE THE INVERSE-DISCRETE-FOURIER-TRANSFORM

y = idft(Y, Fs);
s = idft(S, Fs);

%% PLOT OR PlAY

%time plot
t = (1:length(sig_clean))*T;
figure;
plot(t,sig_clean, 'r');
hold on

% frequency plot
% Yavg = mean(S,2);
% Py = abs(Yavg);
% x = T*(1:length(Py));
% plot(x,Py)

% sound( s (10e3:100e3), Fs)

%% ERRORS

sig_clean = sig_clean(1 : length(s) );
abs_s = mod(sig_clean, 2);
err = (sig_clean - s'); % error due to the STFT
t = (1:length(err))*T;
plot(t, err)
