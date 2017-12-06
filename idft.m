function x = idft(X,Fs)
    % consts
    win_time = 0.020;           %20ms time window
    overlap_ratio = 0.25;       %25% overlap
   
    win_length = 2^nextpow2(win_time*Fs);    %win_length in sample
    overlap = overlap_ratio*win_length;
    nfft = win_length;    
    
    x = istft(X, win_length, overlap, nfft, Fs);
    x = x.'; 
end