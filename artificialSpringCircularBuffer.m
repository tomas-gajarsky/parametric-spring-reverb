%% Artificial Spring with circular buffers 
% - based on Parametric Spring Reverberation Effect
% - the low-frequency chirps block removed
% - the high-chirps block stays the same

clc
clear 

%% User interaction parameters
g_mod = 8; % modulation of delay lines 
g_wet = 0.95; % gain of wet signal
fc = 4300;%the maximum frequency of the first sequence of pulses in Hz
Td = 56; % delay length in miliseconds

%% Parameters
fs = 44100; % sampling frequency in Hz

% low chirps cascaded stretched all pass filters
M = 100; % number of filters in cascade
K = fs/(2 * fc); % stretching factor
K1 = round(K)-1; % rounded stretching factor
d = K - K1; % delay parameter
aa1 = 0.75; % coefficient a1

% delay line modulation
aint = 0.93; % filter coefficient a_int

% high chirps cascaded stretched all pass filters
ah = 0.6; % filter coefficient a_high
Mh = 200;

% general parameters
g_hf = -0.77; % feedback gain of high chirps block

% dry/wet mix
g_dry = 1 - g_wet; % gain of input signal


%% Input noise impluse

xin = rand(1,fs/20); % make a short noise impulse signal
xin = xin - mean(xin); % get rid of DC 
xin = [xin zeros(1,fs*4)]; % pad with zeros

x = xin;

%% Prepare read and write positions, arrays and matrices for blocks
yhf = 0;
y5 = zeros(Mh, 2);
y5l = zeros(1, 2*fs);

y_out = zeros(1, length(x));

x_n = 0;
y_n = zeros(1, 2);

y5_n = 2;
y5_d = 1;

y5l_n = 2;
y5l_d = 1;

y_n_n = 2;
y_n_d = 1;

y1_n = 2;
y1_d = 1;


%% Signal Processing Algorithm of Parametric Spring Reverb (High chirps only)
% low chirp block = C_lf
% high chirp block = C_hf

for n = 1:(length(x))
    
    % white noise generator for delay line modulation
    x_n=rand();
    
    % white noise filtering for delay line modulation
    y_n(y_n_n) = x_n - aint*x_n + aint*y_n(y_n_d);
    
    % set the lengths of multitap delay after modulation  
    % total length plus modulation [samples]
    L = round(Td/1000 * fs - (K*M*((1-aa1)/(1+aa1))) + g_mod*y_n(y_n_n));
    
    % C_hf - delay line
    Lh = round(L/2.3); % length of C_hf delay line

    % update the "read pointer" after delay line modulation
    y5l_dLh = round(mod((y5l_n-Lh+length(y5l)), length(y5l)))+1;
    
    % set the input of the C_hf all-pass cascade matrix 
    y5(1,y5_n) = x(n) - g_hf*yhf;
    
    % C_hf - all-pass cascade
    for k = 2:Mh
       y5(k,y5_n) = ah*y5(k-1,y5_n) + y5(k-1,y5_d) - ah*y5(k,y5_d);
    end  
    y5l(y5l_n) = ah*y5(Mh,y5_n) + y5(Mh,y5_d) - ah*y5l(y5l_d);
    
    % feedback delay line
    yhf = y5l(y5l_dLh);

    % output of the system - block summation
    y_out(n) = g_dry*x(n) + g_wet*y5l(y5l_n);
    
    % increment the "write and read pointers"
    y5_n = y5_n+1;
    y5_d = y5_d+1;
    
    y5l_n = y5l_n+1;
    y5l_d = y5l_d+1;
    
    y_n_n = y_n_n+1;
    y_n_d = y_n_d+1;

    % reset the position of the "write and read pointers"
    if y5_n > 2
       y5_n = 1;
    end
    
    if y5_d > 2
       y5_d = 1;
    end
    
    if y5l_n > length(y5l)
       y5l_n = 1;
    end
    
    if y5l_d > length(y5l)
       y5l_d = 1;
    end
    
    if y_n_n > 2
       y_n_n = 1;
    end
    
    if y_n_d > 2
       y_n_d = 1;
    end
    
end

%% Play the result
soundsc(y_out, fs)

%% Plot the input and output signals
figure
plot(2*x,'r')
hold on
plot(y_out, 'b')
hold off
title('The signal before and after')
xlabel('Time [samples]')
ylabel('Amplitude')
legend('input','output')


%% Write to the audio files
filename_in = 'artificial_spring_in.wav';
filename_out = 'artificial_spring_out.wav';
audiowrite(filename_in,0.5*x,fs);
audiowrite(filename_out,0.5*y_out,fs);
