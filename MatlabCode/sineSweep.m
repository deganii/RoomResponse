function [x] = sineSweep(amplitude, Fs, fjInit, fjFact, J, duration)
%sineSweep: Generation of sine-sweep audio file
%   Generation of sine-sweep audio file

nj = duration * Fs;

for j = 1 : J
    fj = fjInit * fjFact^(j-1);
    x(j,:) = sinWave(amplitude, Fs, fj, duration);
end