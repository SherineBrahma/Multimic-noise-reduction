%% CONSTANTS

% sampling information
Fs = 16000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    

% consts
win_time = 0.020; %sec
overlap_ratio = 0.5;

% lengths
win_len = ceil(win_time.*Fs);
step_len = floor(overlap_ratio.*win_len);

% vectors
win = hann(win_len);

%% LOAD AND STORE THE DATA FILES

S1 = load('data.mat', 'Clean');
S2 = load('data.mat', 'Data');

sc = S1.Clean; % clean signal (1 channel)
sr = S2.Data(:,1); % raw signal (16 channels)

L = size(sc, 1);            % Length of signals (same for Sc and Sr)
t = (0:L-1)*T;              % Time vector

%% PLAY THE AUDIO FILES
% sound(Sr(:,1), Fs);
% disp('Press any button to stop playback...');
% pause;
% clear sound;

%% FFT 
Sc = fft(sc);
Sr = fft(sr);

%stupid edit


if (1) %overall plot
    
    figure(1)
    subplot(2,2,1);
    plot(t(1:L),sc(1:L), 'r')
    title('Clean data')
    xlabel('t (seconds)')
    ylabel('Sc(t)')

    subplot(2,2,2);
    plot(t(1:L), sr(1:L), 'b');
    title('Raw data')
    xlabel('t (seconds)')
    ylabel('Sr(t)')

    P_seg_Sc = abs(Sc/L);
    P_seg_Sc = P_seg_Sc(1:floor(L/2)+1); %only the right side of the FFT is selected
    P_seg_Sc(2:end-1) = 2*P_seg_Sc(2:end-1); %the amplitudes are amplified by a factor 2

    f = Fs*(0:(L/2))/L;

    subplot(2,2,3)
    plot(f,P_seg_Sc,'r') 
    title('Single-Sided Amplitude Spectrum of pc(t)')
    xlabel('f (Hz)')
    ylabel('|PSc(f)|')

    PSr = abs(Sr/L);
    PSr = PSr(1:floor(L/2)+1); %only the right side of the FFT is selected
    PSr(2:end-1) = 2*PSr(2:end-1); %the amplitudes are amplified by a factor 2

    subplot(2,2,4);
    plot(f,PSr, 'b') 
    title('Single-Sided Amplitude Spectrum of pr(t)')
    xlabel('f (Hz)')
    ylabel('|PSr(f)|')
end;

cursor_in = 300*step_len;

%% DISPLAY THE 'BREAK' COMMAND TO EXIT THE FOLLOWING 'WHILE' LOOP
DlgH = figure(2);
H = uicontrol('Style', 'PushButton', ...
                    'String', 'Break', ...
                    'Callback', 'delete(gcbf)');


figure(2);
while (cursor_in < (length(sc)-win_len-step_len)) && ishandle(H) % fragment plot
    %n = win_len %look at 20e-3 window
    ts = cursor_in;
    te = ( ts : ts+win_len-1)';    
    time_seg_sc = sc(cursor_in:(cursor_in+win_len-1)); 
    time_seg_sr = sr(cursor_in:(cursor_in+win_len-1)); 
    
    
    subplot(2,2,1)
    plot(te, time_seg_sc, 'r');
    axis([-inf inf -10 10])
    title('Time domain fragment of sc(t)')
    xlabel('t (s)')
    ylabel('sc(t)')
    
    subplot(2,2,2)
    plot(te, time_seg_sr, 'b');
    axis([-inf inf -10 10])
    title('Time domain fragment of sr(t)')
    xlabel('t (s)')
    ylabel('sr(t)')
                               
    %fft
    f_seg_sc = fft(time_seg_sc.*win);
    
    P_seg_Sc = abs(f_seg_sc/L);
    P_seg_Sc = P_seg_Sc(1:floor(win_len/2)); %only the right side of the FFT is selected
    P_seg_Sc(2:end-1) = 2*P_seg_Sc(2:end-1); %the amplitudes are amplified by a factor 2

    f = (Fs*(0:(win_len/2 -1))/win_len )';
    subplot(2,2,3)
    semilogx(f,P_seg_Sc, 'r')
    axis([-inf inf 0 2e-3])
    title('FFT of windowed sc(t)')
    xlabel('f (Hz)')
    ylabel('Sc(f)')
    
    f_seg_sr = fft(time_seg_sr.*win);
    
    P_seg_Sr = abs(f_seg_sr/L);
    P_seg_Sr = P_seg_Sr(1:floor(win_len/2)); %only the right side of the FFT is selected
    P_seg_Sr(2:end-1) = 2*P_seg_Sr(2:end-1); %the amplitudes are amplified by a factor 2

    f = (Fs*(0:(win_len/2 -1))/win_len )';
    subplot(2,2,4)
    semilogx(f,P_seg_Sr, 'b')
    axis([-inf inf 0 2e-3])
    title('FFT of windowed sr(t)')
    xlabel('f (Hz)')
    ylabel('Sr(f)')
    
    cursor_in = cursor_in+step_len;
    pause(0.2);
end

clc
close all;