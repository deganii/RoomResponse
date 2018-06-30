function [x] = sinWave(amplitude, Fs, fj, duration)
%cosWave: Generator of sine wave

deltaT = 1. / Fs;
nj = duration / deltaT;
x = [0 : 1 : nj-1];
x = amplitude * sin(2. * pi * fj * x * deltaT);