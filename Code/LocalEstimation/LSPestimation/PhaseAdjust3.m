function [dx dy] = PhaseAdjust3(P1_r,P2_r,P1_s,P2_s,fx,fy)

%   PhaseAdjust  Estimation local via ajustement des plans de phase
%   appliquée directement sur les signaux et après avoir calculé les
%   signaux analytiques sur les images entières (used with BM and with BDBM
%   at the 1st iteration)
%   [dx dy] = PhaseAdjust(r,s,fx,fy) computes the 2D estimated local
%   translatioans using phase adjustement
%   Inputs : r,s 2D signals
%            fx,fy frequencies along each direction x and y
%   Outputs : dx,dy estimated local translations

%   Author: Adrian Basarab, Mai 2007.


% Reference
%   [1] A. Basarab, P. Gueth, H. Liebgott, P. Delachartre, Two-dimensional 
%       least-squares estimation for motion tracking in ultrasound elastography,
%       IEEE EMBC, 2007.

% Exemple

% x = 1:30; y = 1:30; [X,Y] = meshgrid(x,y);
% fx = 1/10; fy = 1/10;
% % True shifts
% dx = 0.45; dy = 0.45;
% 
% r = cos(2*pi*fx*X) .* cos(2*pi*fy*Y);
% s = cos(2*pi*fx*(X-dx)) .* cos(2*pi*fy*(Y-dy));
% [r1,r2] = ashahn(r); [s1,s2] = ashahn(s);
% P1_r=atan2(imag(r1),real(r1));
% P2_r=atan2(imag(r2),real(r2));
% P1_s=atan2(imag(s1),real(s1));
% P2_s=atan2(imag(s2),real(s2));
% figure; imagesc(x,y,r); title('r')
% figure; imagesc(x,y,s); title('s')
% 
% [dx_estim dy_estim] = PhaseAdjust3(P1_r,P2_r,P1_s,P2_s,fx,fy)




phi_diff1 = P1_r - P1_s;
phi_diff2 = P2_r - P2_s;

pos = find(abs(phi_diff1)<pi & abs(phi_diff2)<pi);

A = phi_diff1(pos);
B = phi_diff2(pos);

dx = (1/(4*pi*fx))*mean(A-B);
dy = (1/(4*pi*fy))*mean(A+B);


