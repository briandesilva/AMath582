% Amath 582 HW 2
% Brian de Silva
% 1422824

% clear all; close all;

%% Part 1
% load handel
% v = y'/2;
% subplot(2,1,1)
% plot((1:length(v))/Fs,v);
% xlabel('Time [sec]');
% ylabel('Amplitude');
% title('Signal of Interest, v(n)');
% 
% % % Play back audio
% % p8 = audioplayer(v,Fs);
% % playblocking(p8);
% 
% % Identify first and last points of v
% v = v(1:end-1);
% 
% L = (length(v)-1) / Fs;
% t1 = (0:length(v))/Fs;
% t = t1(1:end-1);
% n = length(v);
% k = (2*pi/L)*[0:n/2 - 1, -n/2:-1];
% ks = fftshift(k);
% 
% % Look at frequency components of entire signal
% vt = fft(v);
% subplot(2,1,2)
% plot(ks, abs(fftshift(vt)) / max(abs(vt)));
% xlabel('Frequency (k)')
% ylabel('FFT(v)')
% title('FFT of v unfiltered')
% 
% 
% % Create Gabor filter for spectrogram
% width = 0.1;          % Width of filter
% numsteps = 100;        % Number of time steps to take
% tslide = linspace(0,t(end-1),numsteps);             % Time discretization
% spec = zeros(length(tslide),length(v));             % Preallocate space for spectrogram
% 
% % Create multiple Gabor filters for testing
% % { Gaussian, super Gaussian, Mexican Hat, step function}
% filter = {@(x) exp(-width*(x).^2), @(x) exp(-width*(x).^10),@(x) (1-(x/width).^2).*exp(-((x/width).^2)/2), @(x) (x>-width/2 & x< width/2)};
% 
% 
% figure()
% for j=1:length(tslide)
%     g = filter{4}(t-tslide(j));
% %     g = exp(-width*(t - tslide(j)).^10);     % Gabor filter
%     vg = g.*v;
%     vgt = fft(vg);
%     spec(j,:) = abs(fftshift(vgt));             % Store fft in spectrogram
%     
%     % Make some cool movies of what the Gabor transform is accomplishing
%     subplot(3,1,1), plot(t,v,'k',t,g,'r'), title('Gabor filter and signal'), legend('v','Gabor filter')
%     subplot(3,1,2), plot(t,vg,'k'), title('Gabor filter * signal')
%     subplot(3,1,3), plot(ks, abs(fftshift(vgt))/max(abs(vgt))), title('Gabor transformation of signal')
%     drawnow
%     
% end
% 
% % Plot spectrogram
% figure()
% pcolor(tslide,ks,spec.'), shading interp
% colormap('hot')



%% Part II - Mary had a little lamb

clear all; close all;

%% Piano
% music1.wav plays mhall on piano

tr_piano=16;                     % Record time in seconds
% y=wavread('music1').';    % Read music
[y,Fs] = audioread('music1.wav');
y = y.';
% Fs=length(y)/tr_piano;
figure()
subplot(2,1,1)
plot((1:length(y))/Fs,y);       % Plot data
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow

% Play music
% p8 = audioplayer(y,Fs);
% playblocking(p8);

L = tr_piano;
t1 = linspace(0,16,length(y)+1);
t = t1(1:end-1);
n = length(y);
k = (2*pi/L)*[0:n/2 - 1, -n/2:-1];
ks = fftshift(k);

% Look at frequency components of entire signal
yt = fft(y);
subplot(2,1,2)
plot(ks, abs(fftshift(yt)) / max(abs(yt)));
xlabel('Frequency [Hz]'), ylabel('FFT(y)'), title('FFT of y unfiltered (piano)')

% Create Gabor filter for spectrogram
width = 100;                                                      % Width of filter
numsteps = 100;                                               % Number of time steps to take
tslide = linspace(0,t(end-1),numsteps);            % Time discretization
spec = zeros(length(tslide),length(y));               % Preallocate space for spectrogram

% Create spectrogram using Gabor filter
figure()
for j=1:length(tslide)
    g = exp(-width*(t - tslide(j)).^2);     % Gabor filter
    yg = g.*y;
    ygt = fft(yg);
    spec(j,:) = abs(fftshift(ygt));             % Store fft in spectrogram
    
    % Make some cool movies of what the Gabor transform is accomplishing
    subplot(3,1,1), plot(t,y,'k',t,g,'r'), title('Gabor filter and signal'), legend('y','Gabor filter')
    xlabel('Time [sec]'), ylabel('Amplitude')
    subplot(3,1,2), plot(t,yg,'k'), title('Gabor filter * signal')
    xlabel('Time [sec]'), ylabel('Amplitude')
    subplot(3,1,3), plot(ks, abs(fftshift(ygt))/max(abs(ygt))), title('Gabor transformation of signal')
    xlabel('Frequency [Hz]'), ylabel('Magnitude')
    drawnow
    
end


% Plot relevant portion of spectrogram
figure()
pcolor(tslide,ks,log(spec.'+1)), shading interp
axis([0 15 1400 2200])
colormap('hot'), xlabel('Time [sec]'), ylabel('Frequency [Hz]'), title('Log of Piano Spectrogram (zoomed)')

% Plot full spectrogram
figure()
pcolor(tslide,ks,log(spec.'+1)), shading interp
colormap('hot'), xlabel('Time [sec]'), ylabel('Frequency [Hz]'), title('Log of Piano Spectrogram')

%% Recorder
% music2.wav plays mhall on recorder

figure()
tr_rec=14;                      % Record time in seconds
[y,Fs1]=audioread('music2.wav');
y = y.';
Fs=length(y)/tr_rec;
subplot(2,1,1)
plot((1:length(y))/Fs,y);   % Plot data
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');

% % Play music
% p8 = audioplayer(y,Fs);
% playblocking(p8);

L = tr_rec;
t1 = linspace(0,14,length(y)+1);
t = t1(1:end-1);
n = length(y);
k = (2*pi/L)*[0:n/2 - 1, -n/2:-1];
ks = fftshift(k);

% Look at frequency components of entire signal
yt = fft(y);
subplot(2,1,2)
plot(ks, abs(fftshift(yt)) / max(abs(yt)));
xlabel('Frequency [Hz]'), ylabel('FFT(y)'), title('FFT of y unfiltered (recorder)')

% Create Gabor filter for spectrogram
width = 100;                                                    % Width of filter
numsteps = 100;                                             % Number of time steps to take
tslide = linspace(0,t(end-1),numsteps);          % Time discretization
spec = zeros(length(tslide),length(y));             % Preallocate space for spectrogram

% Create spectrogram using Gabor filter
figure()
for j=1:length(tslide)
    g = exp(-width*(t - tslide(j)).^2);     % Gabor filter
    yg = g.*y;
    ygt = fft(yg);
    spec(j,:) = abs(fftshift(ygt));             % Store fft in spectrogram
    
    % Make some cool movies of what the Gabor transform is accomplishing
    subplot(3,1,1), plot(t,y,'k',t,g,'r'), title('Gabor filter and signal'), legend('y','Gabor filter')
    xlabel('Time [sec]'), ylabel('Amplitude')
    subplot(3,1,2), plot(t,yg,'k'), title('Gabor filter * signal')
    xlabel('Time [sec]'), ylabel('Amplitude')
    subplot(3,1,3), plot(ks, abs(fftshift(ygt))/max(abs(ygt))), title('Gabor transformation of signal')
    xlabel('Frequency [Hz]'), ylabel('Magnitude')
    drawnow
    
end


% Plot relevant portion of spectrogram
figure()
pcolor(tslide,ks,log(spec.'+1)), shading interp
axis([0 14 4000 8000])
colormap('hot'), xlabel('Time [sec]'), ylabel('Frequency [Hz]'), title('Log of Recorder Spectrogram (zoomed)')

% Plot full spectrogram
figure()
pcolor(tslide,ks,log(spec.'+1)), shading interp
colormap('hot'), xlabel('Time [sec]'), ylabel('Frequency [Hz]'), title('Log of Recorder Spectrogram')


