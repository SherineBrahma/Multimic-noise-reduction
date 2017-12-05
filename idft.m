function x = idft(X,Fs)
    % consts
    win_time = 0.020;           %20ms time window
    overlap_ratio = 0.5;        %50% overlap
    win_length = win_time*Fs;    %win_length in sample
    overlap = overlap_ratio*win_length;

    N = length(X);            % Length of the signal

    nfft = 2^nextpow2(win_length);

    % [Y,f,t] = spectrogram(sig_noisy,win_hann,overlap,nfft,Fs);

    x = istft(X, win_length, overlap, nfft, Fs);
end