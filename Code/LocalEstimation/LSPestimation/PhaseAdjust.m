function [dx dy] = PhaseAdjust(r,s,fx,fy)

%   PhaseAdjust  Estimation local via ajustement des plans de phase using the cross correlation function.
%   [dx dy] = PhaseAdjust(r,s,fx,fy) computes the 2D estimated local
%   translatioans using phase adjustement
%   Inputs : r,s 2D signals
%            fx,fy frequencies along each direction x and y
%   Outputs : dx,dy estimated local translations

%   Author: Adrian Basarab, january 2007.


% References
%     [1] A. Basarab, C. Grava, V. Buzuloiu, P. Delachartre, Estimation de
%     décalages subpixéliques par ajustement de la phase des signaux
%     complexes, Gretsi, 2007


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
% [dx_estim dy_estim] = PhaseAdjust(r,s,fx,fy)


lagsX = -(size(r,2)-1):(size(r,2)-1);
lagsY = -(size(r,1)-1):(size(r,1)-1);

R = xcorr2(s,r);

%répérer le maximum de corrélation et sa position
[vect_max vect_pos] = max(R);
[reel_max x_max] = max(vect_max);
y_max = vect_pos(x_max);

[R1,R2]=ashahn(R);

P1=atan2(imag(R1),real(R1));
[dP1odx dP1ody] = gradient(P1);
P2=atan2(imag(R2),real(R2));
[dP2odx dP2ody] = gradient(P2);

[x1 x2 y1 y2] = plane(P1,P2,x_max,y_max);

try
    P1_plane = P1(y1:y2,x1:x2);
    P2_plane = P2(y1:y2,x1:x2);

    x1 = lagsX(x1);
    x2 = lagsX(x2);
    y1 = lagsY(y1);
    y2 = lagsY(y2);
catch
    x1
    x2
    y1
    y2
    figure; mesh(P1)
    figure; mesh(P2)
end

dx = (x1+x2)/2 - (1/(4*pi*fx*(x2-x1+1)*(y2-y1+1)))*sum(sum(P1_plane-P2_plane));
dy = (y1+y2)/2 - (1/(4*pi*fy*(x2-x1+1)*(y2-y1+1)))*sum(sum(P1_plane+P2_plane));
