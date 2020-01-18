function def_grid(im,imbin,pts_ax,pts_lat,pixely,pixelx,dep_ax,dep_lat)

x = pts_lat(1):pts_lat(end);
y = pts_ax(1):pts_ax(end);
[X Y] = meshgrid(x,y);

X_def = X-dep_lat; Y_def = Y-dep_ax;

imbin_def=interp2(X,Y,imbin(floor(y),floor(x)),X_def,Y_def);
figure; imagesc(x*pixelx,y*pixely,im(floor(y),floor(x)).*imbin_def,[0 255]); colormap gray; %axis image
xlabel('Lateral distance [mm]','FontSize',14)
ylabel('Axial distance [mm]','FontSize',14)
figure; imagesc(x*pixelx,y*pixely,imbin(floor(y),floor(x))); colormap(flipud(gray)); axis image
figure; imagesc(x*pixelx,y*pixely,imbin_def); colormap(flipud(gray)); axis image