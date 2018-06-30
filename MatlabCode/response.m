function [ Csin, Ccos ] = response( Y, Fs, fjInit, fjFact, J, duration )
%response: Compute the response to the sine wave

    nj = duration * Fs;

    for j = 1 : J   
        fj = fjInit * fjFact^(j-1);
        Csin(j) = 1. / nj * sum(Y(j,:) .* sinWave(1, Fs, fj, duration));
        Ccos(j) = 1. / nj * sum(Y(j,:) .* cosWave(1, Fs, fj, duration));
        
    end
end

