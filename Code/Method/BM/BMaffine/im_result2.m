function im_result2

% if type_correl2==6
% 
% 
% 
% 
% elseif type_correl2==9



load distribut00.mat
load chpB4.mat

true_lat=readimg('D_x4.img',160,160);
true_ax=readimg('D_y4.img',160,160);
mask=imread('mask.jpg');
mask=rgb2gray(mask);
mask=double(mask);
mask=mask>20;
[y_tail x_tail]=size(true_ax);
[y_tail2 x_tail2]=size(mask);
[x_grid1 y_grid1]=meshgrid(x_tail2/x_tail:x_tail2/x_tail:x_tail2,y_tail2/y_tail:y_tail2/y_tail:y_tail2);
[x_grid2 y_grid2]=meshgrid(1:x_tail2,1:y_tail2);
mask=interp2(x_grid2,y_grid2,mask,x_grid1,y_grid1,'cubic');



% true_lat=readimg('disp_x_0001.img',160,160);
% true_ax=readimg('disp_y_0001.img',160,160);

decim_x = 2;
decim_y = 2;

x=10:150;
y=10:150;

true_lat=true_lat(y,x);
true_ax=true_ax(y,x);
mask=mask(y,x);
pixel_x = 0.1; %0.0752;
pixel_y = 0.025; %0.0196;

x2 = min(scat_x):pixel_x:max(scat_x);
y2 = min(scat_y):pixel_y:max(scat_y);
[X2 Y2] = meshgrid(x2,y2);
[X Y] = meshgrid(x,y);

true_lateral=interp2(X,Y,true_lat,X2,Y2);
true_lateral = true_lateral(1:decim_y:end,1:decim_x:end);

true_axial=interp2(X,Y,true_ax,X2,Y2);
true_axial = true_axial(1:decim_y:end,1:decim_x:end);


mask=interp2(X,Y,mask,X2,Y2);
mask=mask(1:decim_y:end,1:decim_x:end);

a=y_grid(end);
b=x_grid(end);
true_axial = true_axial(y_grid(1):a,x_grid(1):b);
true_lateral = true_lateral(y_grid(1):a,x_grid(1):b);
mask=mask(y_grid(1):a,x_grid(1):b);
mask=mask>0.2;
pixel_x = 0.1*decim_x; 
pixel_y = 0.025*decim_y;
true_lateral=flipud(true_lateral);
true_axial=fliplr(true_axial);
true_lateral(isnan(true_lateral))=0;
true_axial(isnan(true_axial))=0;


axial_error=(true_axial-chp_dep_int(:,:,1).*pixel_y).*mask;
lateral_error=(true_lateral-chp_dep_int(:,:,2).*pixel_x).*mask;



[W1 W2]=size(axial_error);
mask=reshape(mask,W1*W2,1);
ind1=(mask==0);
ax_error=reshape(axial_error,W1*W2,1);
ax_error(ind1)=[];

[W1 W2]=size(lateral_error);
mask=reshape(mask,W1*W2,1);
ind1=(mask==0);
lat_error=reshape(lateral_error,W1*W2,1);
lat_error(ind1)=[];

figure;
subplot(2,1,1)
hist(lat_error,200)
title('histogramme erreur estimAtion(affine)')

subplot(2,1,2)
% hist(ax_error,800)
hist(ax_error,200)
title('histogramme erreur estimAtion(affine)')
% axis([-0.06 0.06 0 10000])

format compact
'moyenne erreur axial estimation avec PBM+model affine'
R(1)=mean(ax_error(:))
'ecart type erreur axial estimation avec PBM+model affine'
R(2)=std(ax_error(:))
'moyenne erreur lateral estimation avec PBM+model affine'
R(3)=mean(lat_error(:))
'ecart type erreur lateral estimation avec PBM+model affine'
R(4)=std(lat_error(:))


save result R
%% affichage

x=pts_lat(1):pts_lat(end);
y=pts_ax(1):pts_ax(end);

figure;subplot(2,2,1);imagesc(x*pixel_x,y*pixel_y,true_lateral);colorbar
xlabel('Lateral distance [mm]','FontSize',14)
ylabel('Axial distance [mm]','Fontsize',14)
title('true lateral motion')
subplot(2,2,2);imagesc(x*pixel_x,y*pixel_y,chp_dep_int(:,:,2).*pixel_x);colorbar
xlabel('Lateral distance [mm]','FontSize',14)
ylabel('Axial distance [mm]','Fontsize',14)
title('lateral estimation PBM+affine')

% title('lateral estimation PBM+translational')

subplot(2,2,3);imagesc(x*pixel_x,y*pixel_y,true_axial);colorbar
xlabel('Lateral distance [mm]','FontSize',14)
ylabel('Axial distance [mm]','Fontsize',14)
title('true axial motion')
subplot(2,2,4);imagesc(x*pixel_x,y*pixel_y,chp_dep_int(:,:,1).*pixel_y);colorbar
xlabel('Lateral distance [mm]','FontSize',14)
ylabel('Axial distance [mm]','Fontsize',14)
title('axial estimation PBM+affine')
% title('axial estimation PBM+tranlational')


figure;subplot(2,1,1);imagesc(x*pixel_x,y*pixel_y,axial_error)
xlabel('Lateral distance [mm]','FontSize',14)
ylabel('Axial distance [mm]','Fontsize',14)
title('axial error estimation')
subplot(2,1,2);imagesc(x*pixel_x,y*pixel_y,lateral_error)
xlabel('Lateral distance [mm]','FontSize',14)
ylabel('Axial distance [mm]','Fontsize',14)
title('lateral error estimation')























