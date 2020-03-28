%% Clear workspace
close all; clear all; clc
%% Load the Handel's Messiah music
load handel
v = y'/2;
t = (1:length(v))/Fs;
L = t(end);
N = length(t);
figure(1);
plot(t, v);
xlabel("Time [sec]");
ylabel("Amplitude");
title("Handel music: Signal of interest, v(t)")

%% Shannon Window (Step-function)
figure(2);
widths = [5 2 1 0.1];
rec=@(x,a,b) ones(1,numel(x)).*(t<(a+b/2) & t>(a-b/2));
for j = 1:length(widths)
    shannon = rec(t, 4, widths(j));
    subplot(4, 1, j)
    plot(t,v,'k'), hold on
    plot(t,shannon, 'r', 'Linewidth', [2])
    set(gca, 'Fontsize', [4])
    ylabel('v(t), g(t)')
end
xlabel('time (t)')

%% Gaussian Window
figure(3);
widths = [0.2 1 5 50];
for j = 1:length(widths)
    gauss = exp(-widths(j)*(t-4).^2);
    subplot(4, 1, j)
    plot(t,gauss, 'r', 'Linewidth', [2])
    set(gca, 'Fontsize', [4])
    ylabel('v(t), g(t)')
end
xlabel('time (t)')

%% Mexican Hat Window
figure(4);
widths = [0.2 1 25 100];
for j = 1:length(widths)
    hat = (1-widths(j)*(t-4).^2).*exp(-(widths(j)*(t-4).^2)/2);
    subplot(4, 1, j)
    plot(t,hat, 'r', 'Linewidth', [2])
    set(gca, 'Fontsize', [4])
    ylabel('v(t), g(t)')
end
xlabel('time (t)')
%% Gabor Transform using Gaussian window
figure(5);
spectrogram = [];
tslide = [0:0.01:L];
k = (2*pi)/(L)*[0:N/2 -N/2:-1];     % Odd N
ks = fftshift(k);
for j = 1:length(tslide)
    gauss = exp(-100*(t-tslide(j)).^2);     % width~0.25sec
    vf = gauss.*v;
    vft = fft(vf);
    spectrogram = [spectrogram; abs(fftshift(vft))];
    subplot(3,1,1);
    plot(t,v,'k',t,gauss,'r',"Linewidth",[2]);
    xlabel("Time [sec]"); ylabel("Amplitude")
    subplot(3,1,2), plot(t,vf,'k')
    xlabel("Time [sec]"); ylabel("Amplitude")
    subplot(3,1,3), plot(ks, abs(fftshift(vft))/max(abs(vft)), 'k')
    xlabel("Frequency [Hz]"); ylabel("Amplitude")
    drawnow
    pause(0.1)
end
%% Gabor Transform using Mexican hat function
figure(6);
spectrogram = [];
tslide = [0:0.05:L];
k = (2*pi)/(L)*[0:N/2 -N/2:-1];     % Odd N
ks = fftshift(k);
for j = 1:length(tslide)
    hat = (1-100*(t-tslide(j)).^2).*exp(-(100*(t-tslide(j)).^2)/2);
    vf = hat.*v;
    vft = fft(vf);
    spectrogram = [spectrogram; abs(fftshift(vft))];
    subplot(3,1,1), plot(t,v,'k',t,hat,'r',"Linewidth",[2])
    subplot(3,1,2), plot(t,vf,'k')
    subplot(3,1,3), plot(ks, abs(fftshift(vft))/max(abs(vft)), 'k')
    drawnow
    pause(0.1)
end

%% Gabor Transform using Shannon function (rectangular function)
figure(7);
spectrogram = [];
tslide = [0:0.05:L];
k = (2*pi)/(L)*[0:N/2 -N/2:-1];      % Odd N
ks = fftshift(k);
for j = 1:length(tslide)
    shannon = rec(t, tslide(j), 0.2);
    vf = shannon.*v;
    vft = fft(vf);
    spectrogram = [spectrogram; abs(fftshift(vft))];
    subplot(3,1,1), plot(t,v,'k',t,shannon,'r',"Linewidth",[2])
    subplot(3,1,2), plot(t,vf,'k')
    subplot(3,1,3), plot(ks, abs(fftshift(vft))/max(abs(vft)), 'k')
    drawnow
    pause(0.1)
end
%% Spectrogram
figure(8);
pcolor(tslide, ks, spectrogram.')
shading interp
colormap(hot)
% axis([2.5 4.5 3000 8000])
xlabel("Time [sec]")
ylabel("Frequency [Hz]")
% title("a = 0.2, FWHM ~ 3.7")
% title("a = 1, FWHM ~ 1.7")
% title("a = 25, FWHM ~ 0.4")
% title("a = 100, FWHM ~ 0.2")


% title("Shannon function")
% title("Mexican Hat function")
title("Gaussian Function")

%% Gabor Transform Plots for paper.
clear all; close all; clc

L=10; n=2048; 
t2=linspace(0,L,n+1); t=t2(1:n); 
k=(2*pi/L)*[0:n/2-1 -n/2:-1]; ks=fftshift(k);
S = sin(10*t).*exp(-(t-5).^2) + ...
    cos(19*t).*tanh(0.5*t).*exp(-(t-2).^2) + ...
    cos(3*t).*tanh(0.5*t).*exp(-(t-8).^2);
St=fft(S);

figure(9); % Time domain 
subplot(2,1,1)
plot(t,S,"k") 
set(gca,"Fontsize",[14]), 
xlabel("Time (t)"), ylabel("S(t)")
title("Non-stationary signal S(t).")
grid on
subplot(2,1,2) % Fourier domain 
plot(ks,abs(fftshift(St))/max(abs(St)),"k"); 
axis([-50 50 0 1]) 
set(gca,"Fontsize",[14]) 
xlabel("frequency (\omega)"), ylabel("FFT(S)")
title("Fourier Transform of the signal.")


