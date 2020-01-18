function tau = UQDE(x1,x2,f1)

%  UQDE  Unbiased Quadrature Delay Estimator
%   tau = UQDE(x1,x2,f1) computes the 2D discrete-time analytic signals 
%   x1, x2 signals, f1 frequency
%   tau is the estimated shift between signals

%   Author: Adrian Basarab, march 2008.

%   References:
%     [1] D.L. Maskell, G. S. Woods, The discret-time quadrature subsample 
%     estimation of delay . IEEE Trans. Instrumentation and Measurement,
%     vol.51, n° 1, 2002.
%
%     [2] H.C. So, A comparative study of two discret-time phase delay 
%     estimators . IEEE Trans. Instrumentation and Measurement,
%     vol.54, n° 6, 2005.


% Fe = 1; f1 = 1/5; sx = 20; d1_r = 0.4; d1_s = 0.6; d = d1_s - d1_r
% N = 20; x = 1:1/Fe:N;
% x1 = cos(2*pi*f1*(x-d1_r)).*exp(-pi*((x-10-d1_r)/sx).^2);
% x2 = cos(2*pi*f1*(x-d1_s)).*exp(-pi*((x-10-d1_s)/sx).^2);
% d_UQDE = UQDE(x1,x2,f1)


delta = 1/(4*f1); N = size(x1,2); x=1:N;
x_Q = x+delta; x1_Q = interp1(x+delta,x1,x+floor(delta)); 
x1_Q=x1_Q(2:(end-floor(delta))); x1=x1((2+floor(delta)):end);
x_Q = x+delta; x2_Q = interp1(x+delta,x2,x+floor(delta)); 
x2_Q=x2_Q(2:(end-floor(delta))); x2=x2((2+floor(delta)):end);

Qm1=x1_Q.*x2; Qm2=x1.*x2; Qm3=x1.*x2_Q; Qm4=x1_Q.*x2_Q;

tau = atan((mean(Qm1)-mean(Qm3))/(mean(Qm2)+mean(Qm4)))/(2*pi*f1);

