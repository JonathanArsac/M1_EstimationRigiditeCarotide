function A = enveloppe2D(xr,m,n)
%ENVELOPPE2D  Discrete-time 2D envelope detection via Hilbert transform.
%   A = HILBERT2(Xr) computes the 2D discrete-time enveloppe detection
%   following the mathematical model cos(2*pi*f1*x)*cos(2*pi*f2*y)*A(x,y)
%   If the input Xr is complex, then only the real part is used: Xr=real(Xr).
%
%   enveloppe2D(Xr,M,N) computes the MxN-point Hilbert transform.  Xr is padded
%   zeros if it has less than MxN points, and truncated if it has more.  
%
%   Example:
%     m = 40; n = 50; fx = 1/10; fy = 2/10; sx = 20; sy = 20;
%     x = 1:n; y = 1:m; 
%     [X,Y] = meshgrid(1:m,1:n); 
%     xr = cos(2*pi*fy*Y).*cos(2*pi*fx*X).*exp(-pi*((X-20)/sx).^2).*exp(-pi*((Y-25)/sy).^2);
%     figure; mesh(enveloppe2D(xr)); colormap(gray); title('Detected 2-D enveloppe')
%     env2D = exp(-pi*((X-20)/sx).^2).*exp(-pi*((Y-25)/sy).^2);
%     figure; mesh(env2D); colormap(gray); title('Imposed 2-D enveloppe')
%
%   See also HILBERT, HILBERT2, ASHAHN.

%   Authors: Adrian Basarab, Philippe Delachartre, january 2007.

%   References:
%     [1] T. Bulow, G. Sommer, Hypercomplex signals - a novel extension of
%     the analytical signal to the multidimensional case, IEEE Transactions
%     on Signal Processing, Vol.49, No.11, 2001.
%
%     [2] S. L. Hahn, Multidimensional complex signals with single-orthant
%     spectra, Proceedings of the IEEE, Vol.80, No.8, August 2002. 

if nargin<2, n=[]; end
if ~isreal(xr)
  warning('enveloppe2D ignores imaginary part of input.')
  xr = real(xr);
end
if isempty(n)
  [m,n] = size(xr);  
end
if m<2 || n<2, 
    x = hilbert(xr); % 1D analytic signal
    return; 
end; 
x = fft2(xr,m,n); % mxn-point 2D FFT.
[u,v] = meshgrid(1:n,1:m);
mo2 = (m + rem(m,2))/2; no2 = (n + rem(n,2))/2;
signu = (u<=no2) - (u==1) - (u>no2) + (u==(no2+1)) - rem(n,2)*(u==no2);
signv = (v<=mo2) - (v==1) - (v>mo2) + (v==(mo2+1)) - rem(m,2)*(v==mo2);

h1 = 1 - j*signu.*signv; % [1] eq.(12) for n=2
as1 = ifft2(x.*h1); % total analytical signal
h2 = -j*signu + signv;% [2] eq.27 et eq.28
as2 = ifft2(x.*h2);
A = sqrt(abs(as1).^2 + abs(as2).^2);

