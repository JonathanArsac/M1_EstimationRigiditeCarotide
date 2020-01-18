function x = hilbert2(xr,m,n)
%HILBERT2  Discrete-time 2D analytic signal via Hilbert transform.
%   X = HILBERT2(Xr) computes the 2D discrete-time analytic signal 
%   X = Xr + i*Xi such that Xi is the Hilbert transform of real image Xr.
%   If the input Xr is complex, then only the real part is used: Xr=real(Xr).
%
%   HILBERT2(Xr,M,N) computes the MxN-point Hilbert transform.  Xr is padded
%   zeros if it has less than MxN points, and truncated if it has more.  
%
%   Example:
%   xr = zeros(40);
%   [X,Y] = meshgrid(1:20,1:20);
%   xr(11:30,11:30) = cos(2*pi*0.1*X).*cos(2*pi*0.1*Y);
%   % 2D envelop detection
%   mesh(abs(hilbert2(xr)))
%
%   See also HILBERT, ASHAHN.

%   Author: Philippe Delachartre, december 2006.

%   References:
%     [1] R. R. Read and S. Treitel, The stabilization of two-dimensional 
%     recursive filters via the discrete Hilbert transform, IEEE Trans.
%     Geo. Sci. Electron., vol. GE-11, pp. 153-160, 207, july 1973.
%
%     [2] F. Peyrin, Y. M. Zhu and R. Goutte, Extension of the notion of 
%     analytic signal for multidimensional signal. Applixation to images, 
%     in Signal Processing III: Theories and Applications, I. T. Young et 
%     al., Eds Amsterdam, The Netherlands: North Holland, 1986, pp. 677-
%     680.

if nargin<2, n=[]; end
if ~isreal(xr)
  warning('HILBERT2 ignores imaginary part of input.')
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
h = 1/2*(signu + signv) + 1;
x = ifft2(x.*h);


