function A = enveloppe1D(xr,f,decimation)
%ENVELOPPE2D  Discrete-time 1D envelope detection via spatial delay and Taylor approx.
%   Example 1:
%     m = 40; fx = 1/12; sx = 20; x = 1:m; 
%     xr = (cos(2*pi*fx*x).*exp(-pi*((x-20)/sx).^2))';
%     A = enveloppe1D(xr,fx,1);
%     figure; plot(xr); hold; plot(A,'r'); legend('signal','envelope')
%   Example 2
%      load imUS.mat; figure; imagesc(log(1+0.1*enveloppe1D(im(:,:,5),0.1,3))); colormap gray
%   See also HILBERT, HILBERT2, ASHAHN, ENVELOPPE2D.

%   Authors: Adrian Basarab, Philippe Delachartre, january 2008.

%   References:
%     [1] M. Schlaikjer, J. P. Bagge, O. M. Soerensen, J. A. Jensen, 
%     Trade off study on different envelope detectors for B-mode imaging,
%     Ultrasonic Symposium, 2003

if ~isreal(xr)
  warning('enveloppe1D ignores imaginary part of input.')
  xr = real(xr);
end

% delay depending on the central frequency
d = floor(1/(4*f));

xr_IQ = xr((d+1):end,:);
xr = xr(1:(end-d),:);

k=1;
for i = 1:decimation:size(xr,1)
    for j = 1:size(xr,2)
        if abs(xr(i,j))<abs(xr_IQ(i,j))
            A(k,j) = abs(xr(i,j))/2 + abs(xr_IQ(i,j));
        else
            A(k,j) = abs(xr(i,j)) + abs(xr_IQ(i,j))/2;
        end
    end
    k=k+1;
end