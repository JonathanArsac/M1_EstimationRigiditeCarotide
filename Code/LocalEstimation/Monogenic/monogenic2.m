
% Compute the monogenic signal by convolution in spatial domain.
% Return amplitude ,phase, orientation, frequency.
% Method described in:
% [1] Non-Rigid Ultrasound Image Registration Based on Intensity and Local Phase Information
% Jonghye Woo - Byung-Woo Hong - Chang-Hong Hu - K. Kirk Shung - C.-C. Jay Kuo - Piotr J. Slomka
% [2] A New Extension of Linear Signal Processing for Estimating Local Properties and Detecting Features
% Michael Felsberg and Gerald Sommer
% [3] M. Felsberg, “Optical Flow Estimation From Monogenic Phase”, First
% International Workshop, IWCM 2004.

% The monogenic signal is here computed with a convolution.


function [amp,phi,theta,freq] = monogenic2(im)

im = im / max([max(max(im)),-min(min(im))]);
    
[M,N] = size(im);

kernelSize = 16;

h1 = zeros(kernelSize,kernelSize);
h2 = zeros(kernelSize,kernelSize);

for u=-kernelSize/2:kernelSize/2-1
    for v=-kernelSize/2:kernelSize/2-1
        if u == 0 && v ==0 
            h1(kernelSize/2+1,kernelSize/2+1) = 0;
            h2(kernelSize/2+1,kernelSize/2+1) = 0;
        else
            h1(u+kernelSize/2+1,v+kernelSize/2+1) = -u / (2*pi*(u^2+v^2)^(3/2));
            h2(u+kernelSize/2+1,v+kernelSize/2+1) = -v / (2*pi*(u^2+v^2)^(3/2));
        end
    end
end

p1 = conv2(im,h1);
p2 = conv2(im,h2);
p1 = p1(kernelSize/2+1:M+kernelSize/2,kernelSize/2+1:N+kernelSize/2);
p2 = p2(kernelSize/2+1:M+kernelSize/2,kernelSize/2+1:N+kernelSize/2);

% amplitude
amp = sqrt(im.*im + p1.*p1 + p2.*p2);

% orientation (d'après [2], formules (6) et (7))
theta = atan(p2./p1);

% phase
phi = angle(im + i*(p1.*cos(theta) + p2.*sin(theta)));

% frequence (d'après [3], formule (13))
diff_p1 = diff(p1,1,1);
diff_p2 = diff(p2,1,2);
diff_im1 = diff(im,1,1);
diff_im2 = diff(im,1,2);
diff_p1 = diff_p1(:,1:N-1);
diff_p2 = diff_p2(1:M-1,:);
diff_im1 = diff_im1(:,1:N-1);
diff_im2 = diff_im2(1:M-1,:);

% freq
diff_phi1 = diff(phi(:,1:end-1),1,1);
diff_phi1 = diff_phi1 - 2*pi*(diff_phi1 > pi) + 2*pi*(diff_phi1 < -pi);
diff_phi2 = diff(phi(1:end-1,:),1,2);
diff_phi2 = diff_phi2 - 2*pi*(diff_phi2 > pi) + 2*pi*(diff_phi2 < -pi);

freq = sqrt(diff_phi1.*diff_phi1 + diff_phi2.*diff_phi2) .* sign(diff_phi1);

% pour que la fréquence ait la meme taille que phi et theta
freq(end+1,:) = freq(end,:);
freq(:,end+1) = freq(:,end);
freq(end,end) = freq(end,end-1);

end