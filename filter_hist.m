function cnts = filter_hist(cnts)
% smooth histogram of uint8 data cnts(1:256)
sig = 10;
cnts = reshape(cnts, 1, 256);
cnts = fft(cnts);
ksq = 1:256;
ksq = min([ (ksq-1).^2; (257-ksq).^2]);
fac = exp(-ksq/(2*sig^2));
cnts = cnts .* fac;
cnts = ifft(cnts);
return