%% Parametric Spring
clc
clear 

%% User interaction parameters
g_mod = 8; % modulation of delay lines 
g_low = 0.97; % gain of wet signal
fc = 4300; % the maximum frequency of the first sequence of pulses in Hz
Td = 56; % delay length in miliseconds
%K = 5.123; % stretching factor

%% Parameters
fs = 44100; % sampling frequency in Hz

% low chirps cascaded stretched all pass filters
M = 100; % number of filters in cascade
%fc = fs/(2 * K); % the maximum frequency of the sequence of pulses in Hz
%fc = 4300; % the maximum frequency of the first sequence of pulses in Hz
K = fs/(2 * fc); % stretching factor
K1 = round(K)-1; % rounded stretching factor
d = K - K1; % delay parameter
aa1 = 0.75; % coefficient a1
aa2 = (1 - d)/(1 + d);  % coefficient a2 determines the fractional delay

% DC filter
fcut = 40; % cut-off frequency of DC high-pass filter in Hz
adc = tan((pi/4) - ((pi*fcut)/fs)); % coeffcient a_dc

% multitap delay line
g_ripple = 0.1; % gain of feedforward path in ripple loop
g_echo = 0.1; % gain of feedforward path in echo loop
N_ripple = 0.5; % the number of ripples
L_ripple = round( 2*K*N_ripple ); % ripple filter length in samples

% stretched IIR (approximates the spectral shape of the chirps)
f_peak = 95; % peak frequency in Hz
B = 130; % bandwith in Hz
K_eq = floor(K); % number of unit delays
R = 1 - (pi*B*K_eq)/(fs); % pole radius
p_cos0 = ((1+R^2)/(2*R)) * cos((2*pi*f_peak*K_eq)/(fs)); % sets pole angle
A0 = 1-R^2; % scaling factor 
aeq1 = -2*R*p_cos0; % filter coefficient a_eq1
aeq2 = R^2; % filter coefficient a_eq2

% coefficients of elliptic low pass filter
b1 = 0.002664752967676;
b2 = -0.013434509058057;
b3 = 0.036183137660622;
b4 = -0.065134704800960;
b5 = 0.089024346522334;
b6 = -0.097829543068792;

a1 = -7.952953597940763;
a2 = 29.748452052047290;
a3 = -68.571433326253030;
a4 = 107.5493966631709;
a5 = -119.7162151903184;
a6 = 95.684751335115460;
a7 = -54.204321680605910;
a8 = 20.833974480614565;
a9 = -4.910863196724604;
a10 = 0.540083712167517;

% delay line modulation
aint = 0.93; % filter coefficient a_int

% high chirps cascaded stretched all pass filters
ah = 0.6; % filter coefficient a_high
Mh = 200;

% general parameters
g_lf = -0.8; % feedback gain of low chirps block
g_hf = -0.77; % feedback gain of high chirps block

c1 = 0.1; % feedback gain of high chirps block into low chirps block
c2 = 0; % feedback gain of low chirps block into high chirps block

% dry/wet mix
g_high = g_low/1000; % gain of high chirps block
g_dry = 1 - g_low - g_high; % gain of input signal


%% Input noise impluse

xin = rand(1,fs/20); % make a short noise impulse signal
xin = xin - mean(xin); % get rid of DC 
xin = [xin zeros(1,fs*4)]; % pad with zeros

x = xin;

%% Prepare arrays and matrices for blocks
y = zeros(1, length(x));
yh = zeros(1, length(x));
y1 = zeros(1, length(y));
y2 = zeros(M+1, length(y));
yf = zeros(1, length(y));
yhf = zeros(1, length(y));
y3 = zeros(1, length(y));
y4 = zeros(1, length(y));
y5 = zeros(Mh+1, length(y));
y_out = zeros(1, length(y));

x_n = zeros(1,length(y));
y_n = zeros(1,length(y));

%% Signal Processing Algorithm of Parametric Spring Reverb
% low chirp block = C_lf
% high chirp block = C_hf

for n = 1:(length(x))
    
    % add the signal from feedback loops to the input
    if n == 1
        y(n) = x(n);
        yh(n) = x(n);
    elseif n>1
        y(n) = x(n) - g_lf*yf(n-1) + c1*y5(Mh+1,n-1);
        yh(n) = x(n) - g_hf*yhf(n-1) + c2*y4(n-1);
    end
   
    % C_lf - DC filter
    if n == 1
        y1(n) = 0.5*y(n) + 0.5*adc*y(n);
    elseif n>1
        y1(n) = 0.5*y(n) - 0.5*y(n-1) + 0.5*adc*y(n) - 0.5*adc*y(n-1) +...
           adc*y1(n-1);
    end
    
    % set the input of the C_lf all-pass cascade matrix 
    y2(1,n) = y1(n);
    
    % C_lf - all-pass cascade
    for j = 2:M+1
        if n == 1
            y2(j,n) = aa1 * y2(j-1,n);
            
        elseif n > 1 && n <= K1
            y2(j,n) = aa1 * y2(j-1,n) + aa1*aa2 * y2(j-1,n-1) -...
                aa2 * y2(j,n-1);
            
        elseif  n > K1 && n <= K1+1
            y2(j,n) = aa1 * y2(j-1,n) + aa1*aa2 * y2(j-1,n-1) +...
                aa2 * y2(j-1,n-K1) - aa2 * y2(j,n-1) -...
                aa1*aa2 * y2(j,n-K1);
        elseif n > K1+1
            y2(j,n) = aa1 * y2(j-1,n) + aa1*aa2 * y2(j-1,n-1) +...
                aa2 * y2(j-1,n-K1) + y2(j-1,n-1-K1) - aa2 * y2(j,n-1) -...
                aa1*aa2 * y2(j,n-K1) - aa1 * y2(j,n-1-K1);
        end
    end
    
    % C_lf - white noise generator for delay line modulation
    x_n(n)=rand();
    
    % C_lf - white noise filtering for delay line modulation
    if n == 1
        y_n(n) = x_n(n) - aint*x_n(n);
        
    elseif n > 1
        y_n(n) = x_n(n) - aint*x_n(n) + aint*y_n(n-1);
    end
    
    % C_lf - set the lengths of multitap delay after modulation
    
    % total length plus modulation [samples]
    L = round(Td/1000 * fs - (K*M*((1-aa1)/(1+aa1))) + g_mod*y_n(n));
    L_echo = round(L/5); % preecho length [samples]
    L0 = L - L_echo - L_ripple; % the remaining delay line length [samples]
    
    % C_lf - multitap delay line
    if n <= L0
        yf(n) = 0;
        
    elseif n > L0 && n <= L_ripple+L0
        yf(n) = g_echo*g_ripple*y2(M+1,n-L0);
        
    elseif n > L_ripple+L0 && n <= L_echo+L0
        yf(n) = g_echo*g_ripple*y2(M+1,n-L0) +...
            g_echo*y2(M+1,n-L_ripple-L0);
        
    elseif n > L_echo+L0 && n <= L_echo+L_ripple+L0     
        yf(n) = g_echo*g_ripple*y2(M+1,n-L0) +...
            g_echo*y2(M+1,n-L_ripple-L0) + ...
            g_ripple*y2(M+1,n-L_echo-L0);
    
    elseif n > L_echo+L_ripple+L0
        yf(n) = g_echo*g_ripple*y2(M+1,n-L0) +...
            g_echo*y2(M+1,n-L_ripple-L0) + ...
            g_ripple*y2(M+1,n-L_echo-L0) + y2(M+1,n-L_echo-L_ripple-L0);  
    end
    
    % C_lf - spectral shaping of chirps with stretched IIR filter
    if n <= K_eq
        y3(n) = (A0/2) * y2(M+1,n);
        
    elseif n > K_eq && n <= 2*K_eq
        y3(n) = (A0/2) * y2(M+1,n) - aeq1*y3(n-K_eq);
        
    elseif n > 2*K_eq
        y3(n) = (A0/2) * (y2(M+1,n) - y2(M+1,n-2*K_eq)) -...
            aeq1*y3(n-K_eq) - aeq2*y3(n-2*K_eq);
    end
    
    % C_lf - low pass filter applied to the ouput of the block
    if n == 1
        y4(n) = b1*y3(n);
        
    elseif n == 2
        y4(n) = b1*y3(n) + b2*y3(n-1) - ...
            a1*y4(n-1);
    
    elseif n == 3
        y4(n) = b1*y3(n) + b2*y3(n-1) + b3*y3(n-2) - ...
            a1*y4(n-1) - a2*y4(n-2);

    elseif n == 4
        y4(n) = b1*y3(n) + b2*y3(n-1) + b3*y3(n-2) + b4*y3(n-3) -...
            a1*y4(n-1) - a2*y4(n-2) - a3*y4(n-3);

    elseif n == 5
        y4(n) = b1*y3(n) + b2*y3(n-1) + b3*y3(n-2) + b4*y3(n-3)...
            + b5*y3(n-4)- ...
            a1*y4(n-1) - a2*y4(n-2) - a3*y4(n-3) - a4*y4(n-4);
    
    elseif n == 6
        y4(n) = b1*y3(n) + b2*y3(n-1) + b3*y3(n-2) + b4*y3(n-3) +...
            b5*y3(n-4) + b6*y3(n-5)- ...
            a1*y4(n-1) - a2*y4(n-2) - a3*y4(n-3) - a4*y4(n-4) - a5*y4(n-5);
    
    elseif n == 7
        y4(n) = b1*y3(n) + b2*y3(n-1) + b3*y3(n-2) + b4*y3(n-3) +...
            b5*y3(n-4) + b6*y3(n-5) + b5*y3(n-6) - ...
            a1*y4(n-1) - a2*y4(n-2) - a3*y4(n-3) - a4*y4(n-4) -...
            a5*y4(n-5) - a6*y4(n-6);
    
    elseif n == 8
        y4(n) = b1*y3(n) + b2*y3(n-1) + b3*y3(n-2) + b4*y3(n-3) +...
            b5*y3(n-4) + b6*y3(n-5) + b5*y3(n-6) + b4*y3(n-7) - ...
        a1*y4(n-1) - a2*y4(n-2) - a3*y4(n-3) - a4*y4(n-4) - a5*y4(n-5) -...
        a6*y4(n-6) - a7*y4(n-7);
    
    elseif n == 9
        y4(n) = b1*y3(n) + b2*y3(n-1) + b3*y3(n-2) + b4*y3(n-3) +...
            b5*y3(n-4) + b6*y3(n-5) + b5*y3(n-6) + b4*y3(n-7) +...
            b3*y3(n-8) - ...
        a1*y4(n-1) - a2*y4(n-2) - a3*y4(n-3) - a4*y4(n-4) - a5*y4(n-5) -...
        a6*y4(n-6) - a7*y4(n-7) - a8*y4(n-8);

    elseif n == 10
        y4(n) = b1*y3(n) + b2*y3(n-1) + b3*y3(n-2) + b4*y3(n-3) +...
            b5*y3(n-4) + b6*y3(n-5) + b5*y3(n-6) + b4*y3(n-7) +...
            b3*y3(n-8) + b2*y3(n-9)- ...
        a1*y4(n-1) - a2*y4(n-2) - a3*y4(n-3) - a4*y4(n-4) - a5*y4(n-5) -...
        a6*y4(n-6) - a7*y4(n-7) - a8*y4(n-8) - a9*y4(n-9);

    elseif n > 10
        y4(n) = b1*y3(n) + b2*y3(n-1) + b3*y3(n-2) + b4*y3(n-3) +...
            b5*y3(n-4) + ...
        b6*y3(n-5) + b5*y3(n-6) + b4*y3(n-7) + b3*y3(n-8) + b2*y3(n-9) +...
        b1*y3(n-10) - ...
        a1*y4(n-1) - a2*y4(n-2) - a3*y4(n-3) - a4*y4(n-4) - a5*y4(n-5) -...
        a6*y4(n-6) - a7*y4(n-7) - a8*y4(n-8) - a9*y4(n-9) - a10*y4(n-10);       
    end
    
    % set the input of the C_hf all-pass cascade matrix 
    y5(1,:) = yh;
    
    % C_hf - all-pass cascade
    for k = 2:Mh+1
        if n == 1
            y5(k,n) = ah*y5(k-1,n);
        elseif n > 1
            y5(k,n) = ah*y5(k-1,n) + y5(k-1,n-1) - ah*y5(k,n-1);
        end
    end
    
    % C_hf - delay line
    Lh = round(L/2.3); % length of C_hf delay line
    if n <= Lh
        yhf(n) = 0;
    elseif n > Lh
        yhf(n) = y5(Mh+1,n-Lh);
    end

    % output of the system - block summation
    y_out(n) = g_dry*x(n) + g_low*y4(n) + g_high*y5(Mh+1,n);
    
end

%% Play the result
soundsc(y_out, fs)

%% Plot the input and output signals
figure
plot(x,'r')
hold on
plot(5*y_out, 'b')
hold off
title('The signal before and after')
xlabel('Time [samples]')
ylabel('Amplitude')
legend('input','output')


%% Write to the audio files
filename_in = 'spring_in.wav';
filename_out = 'spring_out.wav';
audiowrite(filename_in,0.5*x,fs);
audiowrite(filename_out,7*y_out,fs);
