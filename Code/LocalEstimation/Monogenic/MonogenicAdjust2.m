
function [dx dy] = MonogenicAdjust2(im1,im2)

% M. Felsberg, “Optical Flow Estimation From Monogenic Phase”, First
% International Workshop, IWCM 2004.

[a1,phi1,theta1,freq1] = monogenic2(im1);
[a2,phi2,theta2,freq2] = monogenic2(im2);

theta = (theta1+theta2)/2;
freq = (freq1+freq2)/2;

dphi = phi1-phi2;
dphii = dphi.*cos(theta);
dphij = dphi.*sin(theta);

cos2 = cos(theta).*cos(theta).*freq;
sin2 = sin(theta).*sin(theta).*freq;
sincos = cos(theta).*sin(theta).*freq;

scos2 = sum(cos2(find(abs(dphi)<pi)));
ssin2 = sum(sin2(find(abs(dphi)<pi)));
ssincos = sum(sincos(find(abs(dphi)<pi)));

sdphi = [sum(dphii(find(abs(dphi)<pi)));sum(dphij(find(abs(dphi)<pi)))];

M = [scos2,ssincos;ssincos,ssin2];

iM = inv(M);
dep = iM * sdphi;


dx = dep(2);
dy = dep(1);

end
