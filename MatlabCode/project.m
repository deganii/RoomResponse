% ************************************************************************
% Violin System Identification Project
% Step 1: Generation of Sine Sweep Signal
% ************************************************************************

clc; clear all; close all;

%% Sine Sweep Input Generation
amplitude = 0.22;           % Amplitude of sine wav
Fs = 44100;                 % Frequency sampling (Hz)
fjInit = 110;               % Frequency of lowest tone (Hz)
fjFact = 2 ^ (1/12);        % Frequency increment (e.g., tone, semitone)
J = 90;                     % Number of tones
duration = 4;               % Duration of each measurements (s)
fileName = 'track.wav';     % Name of the audio track .wav file

disp(['The experiment duration will be ~', num2str(J * duration), ' s'])

X = sineSweep(amplitude, Fs, fjInit, fjFact, J, duration);

Fjs = zeros(1,J);

% TBD if we want to do this
for j=1:J
    fj = fjInit * fjFact^(j-1);
    wavePeriod = 1. / fj;
    numberOfFullWaves = floor(duration / wavePeriod);
    X(j, (numberOfFullWaves * wavePeriod * Fs) + 1 : duration * Fs) = 0.;
    Fjs(j)=fj;
end
Xmean = mean(X,2);          % should be ~0, if full periods
X = X - Xmean;              % adjust only after  switching to full periods

xTrack = reshape(X', 1, J * Fs * duration);
audiowrite(fileName, xTrack, Fs);

[Csin, Ccos] = response(X, Fs, fjInit, fjFact, J, duration);
inputGain = 2. * sqrt (Csin.^2 + Ccos.^2);      % Numerically unstable
inputPhase = radtodeg (angle (1i * Csin + Ccos) - pi );

figure(); 
subplot(1, 2, 1); semilogy(Fjs,inputGain); title("Input Signal Gain");
subplot(1, 2, 2); plot(Fjs,inputPhase); title("Input Signal Phase"); ylim([-180 0]);

%% Output Data Analysis

forwardTrim = round(Fs / 2);       % Time discarded at beginning of each measurement
backwardTrim = round(Fs / 10);     % Time discarded at end of each measurement

[xTrack, Fs] = audioread('track_old.wav');  %TODO: comment if using newly measured data
[yTrack, FsO] = audioread('measured.wav');

%Lining up signals using xcorr highest peak position
[r, lags] = xcorr(yTrack, xTrack);
[argvalue, argmax] = max(r);
delay = lags(argmax);
delay = delay / Fs;           % delay in s (if breaks here, wrong I/O files used)

trimmedDuration = floor(duration * Fs - forwardTrim - backwardTrim);

yTrack = yTrack(round(delay * Fs) : (round(delay * Fs) + floor(duration * J * Fs) - 1));

Xtrimmed = zeros(J, trimmedDuration);
Ytrimmed = zeros(J, trimmedDuration);

for j=0:J-1
    Xtrimmed(j+1, :) = xTrack(j * Fs * duration + forwardTrim : (j+1) * duration * Fs - backwardTrim - 1)';
    Ytrimmed(j+1, :) = yTrack(j * Fs * duration + forwardTrim : (j+1) * duration * Fs - backwardTrim - 1)';
end
xTrack = reshape(Xtrimmed', 1, J * trimmedDuration);
yTrack = reshape(Ytrimmed', 1, J * trimmedDuration);
audiowrite('measured_trimmed2.wav', yTrack, Fs);

%REMOVE AFTER 
fjInit = fjInit * fjFact;

[CsinMeas, CcosMeas] = response(Ytrimmed, Fs, fjInit, fjFact, J, trimmedDuration / Fs);
gainMeas = 2. * sqrt (CsinMeas.^2 + CcosMeas.^2);
phaseMeas = radtodeg (angle (1i * CsinMeas + CcosMeas) - pi);

% Phase unwrapping
% https://www.ljmu.ac.uk/~/media/files/ljmu/about-us/faculties-and-schools/tae/geri/onedimensionalphaseunwrapping_finalpdf
for j=2:J
    deltaPhase = phaseMeas(j) - phaseMeas(j-1);
    if deltaPhase > 180
        phaseMeas(j:J) = phaseMeas(j:J) - 360;
    elseif deltaPhase < -180
        phaseMeas(j:J) = phaseMeas(j:J) + 360;
    end
end


figure(); 
subplot(1, 2, 1); semilogy(Fjs,gainMeas); title("Output Signal Gain");
subplot(1, 2, 2); plot(Fjs,phaseMeas); title("Output Signal Phase"); 

figure();
subplot(2, 1, 1); plot(gainMeas);
subplot(2, 1, 2); plot(yTrack);

% Creation of inverse filter
n = 50;
%f = [0 0.1 0.5 0.7 0.9 1] % Fjs ./ Fs; % need to normalize (rad / sample)
%f = [0, 0.00249,0.00264,0.00279,0.002966,0.0031,0.003329,0.0035, 1 ];
f = [0 , Fjs ./ Fs * 2, 1]; % (1 / max(Fjs ./ Fs))];
%a = [1 0.25 4.1 0.25 1. 1.]; 
% m = [11.74476, 11.74476,15.2698,31.2721,44.7237,42.226,91.291,168.703, 168.703];
m = 1 ./ gainMeas; % inverse 
%m = m ./ max(m);
m = [0, m, 0];
%b = firpm(n, f, a);


[b,a] = yulewalk(n,f,m);
[h,w] = freqz(b,a,128);
figure();plot(f,m,w/pi,abs(h)); title('Frequency Response of YuleWalk');
% fvtool(b,1)



filteredInput = filter(b, a, xTrack);
%filteredInput = filteredInput ./ max(abs(filteredInput));
filteredInput = filteredInput ./ 50;
audiowrite('filteredX.wav', filteredInput, Fs);
figure();plot(filteredInput);title('Filtered Input')
 
filteredOutput = filter(b, a, yTrack);
%filteredOutput = filteredOutput ./ max(abs(filteredOutput));
filteredOutput = filteredOutput ./ 50;
audiowrite('filteredY.wav', filteredOutput, Fs);
figure();
subplot(2, 1, 1); plot(yTrack);title('UnfFiltered Output')
subplot(2, 1, 2); plot(filteredOutput);title('Filtered Output')

filteredInputResh = reshape(filteredInput', [], 90);
xTrackResh = reshape(xTrack', [], 90);
filteredInputResh = filteredInputResh';
xTrackResh = xTrackResh';

ABtest = [];
for j = 1:90
    ABtest = [ABtest, filteredInputResh(j,:)];
    ABtest = [ABtest, xTrackResh(j,:)];
end
figure(); plot(ABtest)
audiowrite('ABtest.wav', ABtest, Fs);



% Apply to song
[song, Fs] = audioread('monkey_island.mp3');
songMono=song(:,1);
figure(); plot(songMono);
audiowrite('monkey_island.wav', songMono, Fs);

filteredSong = filter(b, a, songMono);
filteredSong = filteredSong ./ (max(abs(filteredSong)));
audiowrite('monkey_island_filtered.wav', filteredSong, Fs);