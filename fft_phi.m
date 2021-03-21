function fft_phi()
% get factors straight to compute gaussian filter via sine FFT's
% diff = 4 in imgaussfilt ===  diff = 64./size(phi').^2

phi = zeros(256,384);
[m,n] = size(phi);

for i = 1:m
    for j = 1:n
        phi(i,j) = exp(-((i-m/2)^2 + (j-n/2)^2)/128);
    end
end

rsq0 = sigma(phi);

diff = 4;
phi1 = imgaussfilt(phi, diff);
rsq1 = sigma(phi1);

fprintf('rsq0,1 = %d %d, diff= %d\n', rsq0, rsq1, diff);

diff2 = 64;   % NB putting box size into k [-0.5, 0.5]
[kx,ky] = meshgrid(ind2k(n), ind2k(m));  % NB flip indices for graphics to x,y
ksq = kx.^2 + ky.^2;

fftphi = dst2(phi);
phi2 = idst2( exp(-ksq*diff2) .* fftphi);

fprintf('rsq2= %d, diff2= %d\n',sigma(phi2), diff2);

return


function kk = ind2k(m)
half = floor(m/2);
kk = [1:half, ((half+1):m)-(m+1)]/m; 

function rsq = sigma(phi)
[m,n] = size(phi);
rsq=0; summ=0;
for i = 1:m
    for j = 1:n
        summ = summ + phi(i,j);
        rsq = rsq + ((i-m/2)^2 + (j-n/2)^2)*phi(i,j);
    end
end
rsq = rsq/summ;