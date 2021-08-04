function new_series = ebisuzaki(series)
    % Written by M. Osman Apr. 2017, based on discussion in "Ebisuzaki, W., 1997: 
    % A Method to Estimate the Statistical Significance of a Correlation When the Data Are Serially Correlated. 
    % J. Climate, 10, 2147?2153, doi: 10.1175/1520-0442(1997)010<2147:AMTETS>2.0.CO;2."
    % INPUT: 
    %   series = some evenly spaced time/space series
    % OUTPUT:
    %   new_series = a pseudorandom variable created by randomly varying the
    %   	phase of each series' Fourier series between 0 and 2pi.  The
    %   	resultant randomly generated series retains the power spectrum of
    %   	the original series (and hence, autocorrelation and variance properties), 
    %       but not the distribution of values.
    series = fft(series);
    
    N  = length(series);
    Nf = floor((N-1)/2);
    phases = rand(Nf,1)*2*pi;
    phases = complex(cos(phases),sin(phases));

    series(2:Nf+1) = series(2:Nf+1).*phases;
    series(end:-1:end-Nf+1) = conj(series(2:Nf+1));

    new_series = real(ifft(series));
end