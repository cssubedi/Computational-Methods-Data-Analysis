%% Clean Workspace
close all; clear all; clc;

%% Load song as vector
% Load the piano recording of the song and create a vector
tr_piano=16; % record time in seconds
y_piano=audioread("music1.wav")'; 
Fs=length(y_piano)/tr_piano;
t_piano=(1:length(y_piano))/Fs;
figure(1);
plot(t_piano,y_piano);
xlabel("Time [sec]"); ylabel("Amplitude");
title("Mary had a little lamb (piano)"); 

% Load the recorder recording of the song and create a vector
tr_rec=14; % record time in seconds
y_recorder=audioread("music2.wav")';
Fs=length(y_recorder)/tr_rec;
t_recorder=(1:length(y_recorder))/Fs;
figure(2);
plot(t_recorder,y_recorder);
xlabel("Time [sec]"); ylabel("Amplitude");
title("Mary had a little lamb (recorder)");

%% Gabor Transform of piano recording using Gaussian window.
figure(3);
L = t_piano(end);
N = length(t_piano);
spectrogram = [];
tslide = [0:0.05:L];
k = (2*pi/L)*[0:N/2-1 -N/2:-1];   % Even N
ks = fftshift(k);
for j = 1:length(tslide)
    gauss = exp(-100*(t_piano-tslide(j)).^2);       % Guassian Window
    yf_piano = gauss.*y_piano;
    yft_piano = fft(yf_piano);
    spectrogram = [spectrogram; abs(fftshift(yft_piano))];
    subplot(3,1,1), plot(t_piano,y_piano,'k',t_piano,gauss,'r',"Linewidth",[2])
    subplot(3,1,2), plot(t_piano,yf_piano,'k')
    subplot(3,1,3), plot(ks, abs(fftshift(yft_piano))/max(abs(yft_piano)), 'k')
    drawnow
    pause(0.1)
end

% Logarithmic Spectrogram
Spectrogram
figure(4);
pcolor(tslide, ks, log(spectrogram.' + 1))
axis([0 15 1400 2200])
shading interp
colormap(hot)
xlabel("Time [Sec]")
ylabel("Frequency [Hz]")

% Zoomed Logarithmic Spectrogram
figure(5);
pcolor(tslide, ks, log(spectrogram.' + 1))
shading interp
colormap(hot)
xlabel("Time [Sec]")
ylabel("Frequency [Hz]")

%% Gabor Transform of recorder recording using Gaussian window.
figure(6);
L = t_recorder(end);
N = length(t_recorder);
spectrogram = [];
tslide = [0:0.05:L];
k = (2*pi/L)*[0:N/2-1 -N/2:-1];   % Even N
ks = fftshift(k);
for j = 1:length(tslide)
    gauss = exp(-100*(t_recorder-tslide(j)).^2);    % Gaussian Window
    yf_recorder = gauss.*y_recorder;
    yft_recorder = fft(yf_recorder);
    spectrogram = [spectrogram; abs(fftshift(yft_recorder))];
    subplot(3,1,1), plot(t_recorder,y_recorder,'k',t_recorder,gauss,'r',"Linewidth",[2])
    subplot(3,1,2), plot(t_recorder,yf_recorder,'k')
    subplot(3,1,3), plot(ks, abs(fftshift(yft_recorder))/max(abs(yft_recorder)), 'k')
    drawnow
    pause(0.1)
end

% Logarithmic Spectrogram
figure(7);
pcolor(tslide, ks, log(spectrogram.' + 1))
axis([0 14 4000 8000])
shading interp
colormap(hot)
xlabel("Time [Sec]")
ylabel("Frequency [Hz]")
title("Spectrogram of Recorder recording.")

% Zoomed Logarithmic Spectrogram
figure(8);
pcolor(tslide, ks, log(spectrogram.' + 1))
shading interp
colormap(hot)
xlabel("Time [Sec]")
ylabel("Frequency [Hz]")