function [dx dy] = PhaseAdjust2(r,s,fx,fy,US_IRM)

%   PhaseAdjust  Estimation local via ajustement des plans de phase
%   appliqu�e directement sur les signaux (used only with BDBM after the 1st iteration)
%   [dx dy] = PhaseAdjust(r,s,fx,fy) computes the 2D estimated local
%   translatioans using phase adjustement
%   Inputs : r,s 2D signals
%            fx,fy frequencies along each direction x and y
%            US_IRM we rescale or not the signals
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
% figure; imagesc(x,y,r); title('r')
% figure; imagesc(x,y,s); title('s')
% 
% [dx_estim dy_estim] = PhaseAdjust2(r,s,fx,fy,1)

if US_IRM==0
    r = r - (min(min(r))+max(max(r)))/2; 
    s = s - (min(min(s))+max(max(s)))/2;
end

[R1,R2]=ashahn(r);
P1_r=atan2(imag(R1),real(R1));
P2_r=atan2(imag(R2),real(R2));

[S1,S2]=ashahn(s);
P1_s=atan2(imag(S1),real(S1));
P2_s=atan2(imag(S2),real(S2));

phi_diff1 = P1_r - P1_s;
phi_diff2 = P2_r - P2_s;

pos = find(abs(phi_diff1)<pi & abs(phi_diff2)<pi);

A = phi_diff1(pos);
B = phi_diff2(pos);

dx = (1/(4*pi*fx))*mean(A-B);
dy = (1/(4*pi*fy))*mean(A+B);


