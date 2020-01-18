function [g,d]=compression_orientation(depl_lat,pixel_ax,pixel_lat)

[mini x_min] = min(abs(depl_lat'));
P = polyfit(1:size(x_min,2),x_min,1);
% figure; imagesc(depl_lat); colorbar
% hold; plot(x_min,1:size(x_min,2),'k'); plot(P(1)*(1:size(x_min,2))+P(2),1:size(x_min,2));      
P_mm = polyfit((1:size(x_min,2))*pixel_ax,x_min*pixel_lat,1);
g = P_mm(1);
d = P_mm(2);
