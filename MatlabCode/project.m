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

% TBD if we want to do this
for j=1:90
    fj = fjInit * fjFact^(j-1);
    wavePeriod = 1. / fj;
    numberOfFullWaves = floor(duration / wavePeriod);
    X(j, (numberOfFullWaves * wavePeriod * Fs) + 1 : duration * Fs) = 0.;
end
Xmean = mean(X,2);          % should be ~0, if full periods
X = X - Xmean;              % adjust only after  switching to full periods

xTrack = reshape(X', 1, J * Fs * duration);
audiowrite(fileName, xTrack, Fs);

[Csin, Ccos] = response(X, Fs, fjInit, fjFact, J, duration);
inputGain = 2. * sqrt (Csin.^2 + Ccos.^2);      % Numerically unstable
inputPhase = radtodeg (angle (1i * Csin + Ccos) - pi / 2);

figure(); 
subplot(1, 2, 1); plot(inputGain); title("Input Signal Gain"); xlim([0 J]);
subplot(1, 2, 2); plot(inputPhase); title("Input Signal Phase"); xlim([0 J]);

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

yTrack = yTrack(round(delay * Fs) : (round(delay * Fs) + floor(duration * J * Fs) - 1));

Xtrimmed = zeros(J, floor(duration * Fs - forwardTrim - backwardTrim));
Ytrimmed = zeros(J, floor(duration * Fs - forwardTrim - backwardTrim));

for j=0:J-1
    Xtrimmed(j+1, :) = xTrack(j * Fs * duration + forwardTrim : (j+1) * duration * Fs - backwardTrim - 1)';
    Ytrimmed(j+1, :) = yTrack(j * Fs * duration + forwardTrim : (j+1) * duration * Fs - backwardTrim - 1)';
end
