function b=dst2(arg1,mrows,ncols)
%dst2 2-D discrete sine transform.
%   B = dst2(A) returns the discrete sine transform of A.
%   The matrix B is the same size as A and contains the
%   discrete sine transform coefficients.
%
%   B = dst2(A,[M N]) or B = dst2(A,M,N) pads the matrix A with
%   zeros to size M-by-N before transforming. If M or N is
%   smaller than the corresponding dimension of A, dst2 truncates
%   A. 
%
%   This transform can be inverted using Idst2.
%
%   Class Support
%   -------------
%   A can be numeric or logical. The returned matrix B is of 
%   class double.
%
%   Example
%   -------
%       I = imread('cameraman.tif');
%       J = dst2(I);
%       imshow(log(abs(J)),[]), colormap(jet), colorbar
%
%   The commands below set values less than magnitude 10 in the
%   dst matrix to zero, then reconstruct the image using the
%   inverse dst function Idst2.
%
%       J(abs(J)<10) = 0;
%       K = idst2(J);
%       figure, imshow(I)
%       figure, imshow(K,[0 255])
%
%   See also FFT2, Idst2, IFFT2.

%   Autor: Said BOUREZG  
%   Electronics Engineer   Option:Communication .
%   Date: 06.02.2015
% Adress:                          Said BOUREZG
%                               Elbassatine street
%                                 28038 Tarmount
%                               M'sila --- Algeria 
% Email:  said.bourezg@yahoo.fr
% Mobile: +213 796 018049 
% If you can improve this code furtherly or add arithmatic coding stage,
% please let me know. Thanks
%   Filename: dst2.m (Matlab)
%   Copyright 2015 Said BOUREZG.

[m, n] = size(arg1);
% Basic algorithm.
if (nargin == 1),
  if (m > 1) && (n > 1),
    b = dst(dst(arg1).').';
    return;
  else
    mrows = m;
    ncols = n;
  end
end

% Padding for vector input.
a = arg1;
if nargin==2, ncols = mrows(2); mrows = mrows(1); end
mpad = mrows; npad = ncols;
if m == 1 && mpad > m, a(2, 1) = 0; m = 2; end
if n == 1 && npad > n, a(1, 2) = 0; n = 2; end
if m == 1, mpad = npad; npad = 1; end   % For row vector.

% Transform.

b = dst(a, mpad);
if m > 1 && n > 1, b = dst(b.', npad).'; end
