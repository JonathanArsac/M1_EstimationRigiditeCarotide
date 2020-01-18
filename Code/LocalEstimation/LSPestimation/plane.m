function [x1 x2 y1 y2] = plane(P1,P2,x0,y0)

%PLANE  Plane extraction from txo 2-D analytical signals phases
%   OUTPUT : the limits in x and y of the planes
%   INPUT : phases P1 and P2 and central point of the plane phases to
%   extract (coordinates of correlation function maximum)
%   Planes will be symetrical compared to y0 in vertical direction
%   See also ASHAHN.

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
% R = xcorr2(s,r);
% 
% [vect_max vect_pos] = max(R); [reel_max x_max] = max(vect_max); y_max = vect_pos(x_max);
% [R1,R2]=ashahn(R);
% P1=atan2(imag(R1),real(R1)); P2=atan2(imag(R2),real(R2));
% [x1 x2 y1 y2] = plane(P1,P2,x_max,y_max);
% figure; subplot(2,2,1); mesh(x,y,P1); title('P1')
% subplot(2,2,2); mesh(x,y,P2); title('P2')
% subplot(2,2,3); mesh(x,y,P1(y1:y2,x1:x2)); title('Plane extracted from P1')
% subplot(2,2,4); mesh(x,y,P2(y1:y2,x1:x2)); title('Plane extracted from P2')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 2;
h=3;
l=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dP1odx dP1ody] = gradient(P1);

dP1odx(find(abs(dP1odx)<s)) = l;
dP1odx(find(abs(dP1odx)>s)) = h;

dP1ody(find(abs(dP1ody)<s)) = l;
dP1ody(find(abs(dP1ody)>s)) = h;

[dP2odx dP2ody] = gradient(P2);

dP2odx(find(abs(dP2odx)<s)) = l;
dP2odx(find(abs(dP2odx)>s)) = h;

dP2ody(find(abs(dP2ody)<s)) = l;
dP2ody(find(abs(dP2ody)>s)) = h;

dP12xy = dP1odx.*dP2odx.*dP1ody.*dP2ody;
dP12xy(find(dP12xy>l)) = h;
% figure; imagesc(dP12xy); colormap gray

y1 = y0;
y2 = y0;

aire(1) = 0; aire(2)=0;
j=2;
x1(1) = 0;
x2(1) = inf;

while aire(j)>=aire(j-1)

    x1(j) = x0;
    x2(j) = x0;

    while x1(j)>2 && dP12xy(y1,x1(j))==l && dP12xy(y2,x1(j))==l

        x1(j) = x1(j) - 1;

    end

    while x2(j)<(size(dP12xy,2)-2) && dP12xy(y1,x2(j))==l && dP12xy(y2,x2(j))==l

        x2(j) = x2(j) + 1;

    end

    aire(j+1) = (x2(j)-x1(j))*(y2-y1);

    if x1(j) > x1(1) && aire(j+1)>=aire(j) && x1(j) < x2(1)
        x1(1)=x1(j);
    else
    end

    if x2(j) < x2(1) && aire(j+1)>=aire(j) && x2(j) > x1(1)
        x2(1)=x2(j);
    else
    end

    y1 = y1-1;
    y2 = y2+1;
    j=j+1;
end

y1 = y1+2;
y2 = y2-2;
x1 = x1(1);
x2 = x2(1);
