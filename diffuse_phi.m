function phi = diffuse_phi(phi, colony, source, lam, niter)
% Use sine FFT and gaussian filter in a consistent manner to compute finite time step of 
%   d phi/dt = D lapl phi - lam*phi + source
% where source is time indep. After each step, phi outsize of colony is zeroed to impose
% boundary conditions.  Test code fft_phi() checks equivalence between units in FFT and
% units in imgasssfilt() from matlab, which presumably is faster than my FFT's
%   Alternative coding is simple euler step finite difference, and is too slow.
% For lam=0.01, sig=4, niter = 1000 gets good convergence. For lam=0, sig=4 niter=10000
% converged.

sig_flt = 4;  % gaussian filter radius in pixels after resizing by 2x

source = filt_source(source, lam, 2*sig_flt);
phi = imresize(phi, 1/2);
colony = imresize(colony, 1/2);
source = imresize(source, 1/2);

% phi2 = zeros(size(phi));
% [m,n] = size(phi);
for i = 1:niter
    phi = zero_noncolony(phi, colony);
    phi = exp(-lam)*imgaussfilt(phi, sig_flt) + source;
%     phi2(2:(m-1),:) = (phi(1:(m-2),:) + phi(3:m,:))/4;
%     phi2(:,2:(n-1)) = phi2(:,2:(n-1)) + (phi(:,1:(n-2)) + phi(:,3:n))/4;
%     % split step , lam ~>1 for stability DOES NOT WORK for lam large. 
%     if lam < 0.25
%         phi2 = phi2 - lam*phi + source;
%         phi = max(phi2, 0);
%     else
%         phi = exp(-lam)*phi2 + (1/lam)*(1 - exp(-lam))*source;
%     end
end
phi = imresize(phi, 2);

function phi = zero_noncolony(phi, colony)
% assume adsorbing boundaries around image and colony edge
[m,n] = size(phi);
phi(1,:) = 0;
phi(m,:) = 0;
phi(:,1) = 0;
phi(:,n) = 0;
phi(~colony) = 0;

function phi = filt_source(phi, lam, sig)
% input phi, decay rate, and sigma for gaussian filter
% From fft_phi: sig = 4 in imgaussfilt <===>  diff = 64./size(phi').^2

diff = 4*sig^2;   % NB putting box size into kx,y [-0.5, 0.5]
[m,n] = size(phi);
[kx,ky] = meshgrid(ind2k(n), ind2k(m));  % NB flip indices for graphics to x,y
ksq = diff*(kx.^2 + ky.^2);
fac = ksq + lam;
fac = (1 - exp(-fac))./(fac + 1.e-10);  % solves initial value problem

fftphi = dst2(phi);
phi = idst2(fac .* fftphi);
return

function kk = ind2k(m)
half = floor(m/2);
kk = [1:half, ((half+1):m)-(m+1)]/m; 