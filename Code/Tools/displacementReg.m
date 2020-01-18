function [depl_lat_cumule,depl_ax_cumule] = displacementReg(pts_ax,pts_lat,depl_ax,depl_lat,depl_lat_cumule,depl_ax_cumule)

k = size(depl_lat_cumule,3);

x = pts_lat(1):pts_lat(end);
y = pts_ax(1):pts_ax(end);
[X Y] = meshgrid(x,y);

X_def = X+sum(depl_lat_cumule(:,:,:),3); Y_def = Y+sum(depl_ax_cumule(:,:,:),3);

temp = interp2(X,Y,depl_ax,X_def,Y_def);
temp(find(isnan(temp))) = 0; depl_ax_cumule(:,:,k+1) = temp;
temp = interp2(X,Y,depl_lat,X_def,Y_def);
temp(find(isnan(temp))) = 0; depl_lat_cumule(:,:,k+1) = temp;
