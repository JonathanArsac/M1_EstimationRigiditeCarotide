function [x1,x2] = ashahn(xr,m,n)
%ASHAHN  Discrete-time 2D analytic signal via Hilbert transform.
%   X = ASHAHN(Xr) computes the 2D discrete-time analytic signals 
%   x1 and x2 with xr = 0.5.Real{x1 + X2} and the fourier transforms: 
%   F{x1} = F{xr}(1 + sign(u))(1 + sign(v))
%   F{x2} = F{xr}(1 - sign(u))(1 + sign(v)).
%   If the input Xr is complex, then only the real part is used: Xr=real(Xr).
%
%   ASHAHN(Xr,M,N) computes the MxN-point Hilbert transform.  Xr is padded
%   zeros if it has less than MxN points, and truncated if it has more.  
%
%
%   See also HILBERT, HILBERT2.

%   Author: Adrian Basarab, Philippe Delachartre, december 2006.

%   References:
%     [1] S. L. Hahn, Hilbert transforms in signal processing. Boston, MA:  
%     Artech House, 1996.  
%
%     [2] T. Bülow, G. Sommer, Hypercomplex Signals - A novel extension of
%     the analytic signal to the multidimensional Case. IEEE Trans. Signal
%     Processing, vol. 49, n° 11, 2001, pp. 2844-2852.

if nargin<2, m=[]; n=[]; end
if nargout<2, x1 = hilbert2(xr,m,n); return; end
if ~isreal(xr)
  warning('ASHAHN ignores imaginary part of input.')
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
h1 = (1 + signu).*(1 + signv);  % [2] eq. (16)
x1 = ifft2(x.*h1);
h2 = (1 - signu).*(1 + signv);  % [2] eq. (18)
x2 = ifft2(x.*h2);

