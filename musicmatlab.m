close all
clear
clc

t = [0:.000125:.20];

Al = sin(2*pi*220*t); % The 'l' attached to A is to show it's a note below middle C.
A = sin(2*pi*440*t);
C = sin(2*pi*523.25*t);
D = sin(2*pi*587.33*t);
E = sin(2*pi*659.26*t);
midc = sin(2*pi*261.6*t);
G = sin(2*pi*783.99*t);
Gl = sin(2*pi*196*t); % The 'l' attached to 'G' is to show it is below middle C
B = sin(2*pi*493.88*t);
F = sin(2*pi*349.23*t);
Fl = sin(2*pi*174.61*t); % The 'l' attached to 'F' is to show it is below middle C
Rs = sin(2*pi*000*t); % I use this for a rest note
Dl = sin(2*pi*293.66*t); % The 'l' attached to the 'D' is to show it is below middle c
Gs = sin(2*pi*392*t); %This G note was not below middle c but not a higher not either.
Bl = sin(2*pi*246.94*t); % The 'l' attached to the 'B' is to show it is below middle c
El = sin(2*pi*164.81*t); % Same condition as above.
As = sin(2*pi*466.16*t); % The 's' stands for a # note

% These are notes that needs to be played at the same time.
% Using 'Var' = ('x' + 'y') I was able to do this.
AlA = (Al + A);
DAl = (D + Al);
Cmidc = (C + midc);
GGl = (G + Gl);
CGl = (C + Gl);
FFl = (F + Fl);
BFl = (B + Fl);
% This is where I create my 'measures' and lines.
line1 = [AlA,A,C,A,D,A,E,DAl,Cmidc,C,E,C,G,C,E,Cmidc,GGl,G,B,G,C,G,D,CGl];
line2 = [FFl,F,A,F,B,F,C,BFl];
line3 = [FFl,F,A,F,B,F,C,BFl];

% This is where I create my actual song. I organize my lines in the correct
% order according to the sheet music I used for the actual song.
song  =  [line1,line2,line1,line3];

audiowrite('MortalKombat.wav',song, 9000);
[x, fs] = audioread('MortalKombat.wav');
noise= awgn(x,5);
audiowrite('MortalKombatWithNoise.wav',noise, 9000);
N = length(x);
X = fft(x);
X = X/N;
X1 = fftshift(abs(X));
f = fs*linspace(0,1, N) - fs/2;

figure
% subplot(5, 1, 1)
plot(f,X1);
ylabel('|X(f)|'); xlabel('Frequency (Hz)')
title('Frequency Content of Input Signal')
xlim([-600 600])

figure()
% subplot(5, 1, 1)
plot(f, noise);
xlim([-600 600])
ylabel('noisy signal'); xlabel('Frequency (Hz)')
title('Noisy signal frequency content')

% build filter
% The code below shows an example on how to build and use a digital filter.
cutoff = 700; % cutoff frequency (Hz) for the filter
% build a 5th-order butterworth low pass digital filter
% to keep the music but filter out the noise
% b: numerator coeff of system transfer function, a: denominator coeff
% system transfer function is discussed in section 9.7 of the textbook.
% In the command below, 5 is the order of the filter, the second parameter
% specify the cutoff frequency (with value between 0 and 1), which is
% normalized to half the sample rate, 'low' denotes the type of the filter,
%  low pass filter
[b, a] = butter(5, 2*cutoff/(fs), 'low');
H = freqz(b, a, f, fs); % get frequency response of filter
figure()
plot(f, abs(H)) % plot frequency response (magnitude) of filter
ylabel('|H(f)|'); xlabel('Frequency (Hz)')
title('Filter  Characteristics')

% Signal processing: pass the input through the lowpass filter
y = filter(b, a, noise); % using difference equation
audiowrite('MortalKombatFiltered.wav',y,9000);

% plot frequency content of filter output
N = length(y);
Y = fft(y(1:N))/N; % generate FFT and divide by num of points to normalize FFT;
f = fs*linspace(0,1, N) - fs/2; % generate frequency vector
Y1 = fftshift(abs(Y)); % magnitude of the double-sided spectrum
figure
plot(f,Y1);
ylabel('|Y(f)|'); xlabel('Frequency (Hz)')
title('Frequency Content of Filtered Output Signal')

tt = (1:length(x))/fs;         % time
dur = find(tt>0.1 & tt<0.12);   % set time duration for waveform plot
figure
plot(tt(dur),x(dur))
axis tight
title('Waveform of song in time domain' )

p = 2*abs( Y); % compute power at each frequency
figure
semilogy(f,p)
axis([0 4000 10^-4 1])
title('Power Spectrum of the song ')  % generate spectogram
