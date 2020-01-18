% Compute the difference of gaussion (DoG) of an image.
% Return a bandpass filtered image.
% Parameters:   - im: the input image
%               - sigma1: standard deviation of the 1st filter
%               - sigma2: standard deviation of the 2nd filter
%               - sizeY: size of the filter
%               - sizeX: size of the filter

function imb=dog(im,sigma1,sigma2,sizeY,sizeX)

imb = abs(filter2(fspecial('gaussian',[sizeY,sizeX],sigma1),im) - ...
            filter2(fspecial('gaussian',[sizeY,sizeX],sigma2),im));
imb = imb - mean2(imb)*ones(size(imb));

end