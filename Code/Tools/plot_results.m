function plot_results(image1,axial_motion,lateral_motion,axial_coord,lateral_coord,pixel_size_ax,pixel_size_lat)

% Plot the results of motion estimation (axial and lateral motion maps,
% motion vectors superimposed to the reference image, axial strain image)

%plot_results(image1,axial_motion,lateral_motion,axial_coord,lateral_coord,pixel_size_ax,pixel_size_lat)

%INPUT
%  image1 : Image1 (B-mode)
%  axial_motion : axial component of estimated motion dense field
%  lateral_motion : lateral component of estimated motion dense field
%  lateral_coord  :  lateral positions of estimated pixels
%  axial_coord  :  axial positions of estimated pixels
%  pixel_size_ax : axial pixel size in mm
%  pixel_size_lat : lateral pixel size in mm

% Adrian Basarab, 2007



siz = 14; % police size for axis labels
pas_x = 1; % lateral step of plotting the motion vectors
pas_y = 80; % axial step of plotting the motion vectors


axial_motion = axial_motion*pixel_size_ax;
lateral_motion = lateral_motion*pixel_size_lat;

x=lateral_coord(1):lateral_coord(end);
y=axial_coord(1):axial_coord(end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimated zone (red rectangle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;imagesc((1:size(image1,2))*pixel_size_lat,(1:size(image1,1))*pixel_size_ax,image1);colormap(gray);
rectangle('Position',[x(1)*pixel_size_lat,y(1)*pixel_size_ax,(x(end)-x(1))*pixel_size_lat,(y(end)-y(1))*pixel_size_ax],'edgecolor','r')
xlabel('Lateral distance [mm]','FontSize',14)
ylabel('Axial distance [mm]','Fontsize',14)
%axis image; 
title('Image 1 and red rectangle corresponding to processed region')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Motion axial and lateral maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; imagesc(x*pixel_size_lat,y*pixel_size_ax,axial_motion); colorbar; mp = ColorSpiral; colormap(mp); title('Axial motion map [mm]'); 
xlabel('Lateral distance [mm]','FontSize',siz)
ylabel('Axial distance [mm]','Fontsize',siz)
% axis image
figure; imagesc(x*pixel_size_lat,y*pixel_size_ax,lateral_motion); colorbar; mp = ColorSpiral; colormap(mp); title('Lateral motion map [mm]'); 
% axis image
xlabel('Lateral distance [mm]','FontSize',siz)
ylabel('Axial distance [mm]','Fontsize',siz)
% axis image


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Motion vectors superimposed to the reference image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depl_ax_quiver = axial_motion(1:pas_y:end,1:pas_x:end);
depl_lat_quiver = lateral_motion(1:pas_y:end,1:pas_x:end);

[ax_x_quiver ax_y_quiver] = meshgrid(floor(x),floor(y));


image_init_env = image1(floor(y),floor(x));
figure; imagesc(floor(x)*pixel_size_lat,floor(y)*pixel_size_ax,image_init_env); colormap(gray); 
hold on
quiver(floor(ax_x_quiver(1:pas_y:end,1:pas_x:end))*pixel_size_lat,floor(ax_y_quiver(1:pas_y:end,1:pas_x:end))*pixel_size_ax,depl_lat_quiver,depl_ax_quiver,1,'g');colormap(gray);
xlabel('Lateral distance [mm]','FontSize',siz)
ylabel('Axial distance [mm]','Fontsize',siz)
title('Motion vectors superimposed to the reference image')
% axis image


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Axial strain image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nc = floor(0.15*[size(axial_motion,1)])+1; % habituellement la longueur du filtre dérivateur est 15% de la taille de l'image dans la direction de dérivation
eyy = derivate2(axial_motion/pixel_size_ax,Nc)*100;
figure; imagesc(x*pixel_size_lat,((y(1)+floor(Nc/2)):((y(end)-floor(Nc/2))))*pixel_size_ax,eyy); colorbar; colormap(flipud(gray)); mp = ColorSpiral; colormap(mp);
xlabel('Lateral distance [mm]','FontSize',siz)
ylabel('Axial distance [mm]','Fontsize',siz)
title('Axial strain image [%]')
% axis image
