clear; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill in the path of images location and the result file name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%schema = imread('parameters.JPG'); figure; imagesc(schema); axis image


% Name of the result file
resultat = 'result';

decim_x = 2;
decim_y = 1;

load adrian_0_sta2_result.mat
im(:,:,1) = image_finale(1:decim_y:end,1:decim_x:end);

load adrian_1_sta2_result.mat
im(:,:,2) = image_finale(1:decim_y:end,1:decim_x:end);


%% IMAGE PARAMETERS %%

% normalized frequency vector
%     fy   fx
f = [0.12 0.1].*[decim_y decim_x];

% pixel size in mm
pixel_size_ax = 0.0196*decim_y;
pixel_size_lat = 0.15*decim_x; %Pour l'echographe Ultrasonix taille pixel 40 mm (taille sonde) /128 lignes = 0.3125 mm



%%% images log B-scan %%%
log_im(:,:,1) = log(abs(hilbert(im(:,:,1))));
log_im(:,:,2) = log(abs(hilbert(im(:,:,2))));

figure;imagesc((1:size(im,2))*pixel_size_lat,(1:size(im,1))*pixel_size_ax,log_im(:,:,1));colormap(gray);
xlabel('Lateral distance [mm]','FontSize',14)
ylabel('Axial distance [mm]','Fontsize',14)
axis image; title('Image 1')
figure;imagesc((1:size(im,2))*pixel_size_lat,(1:size(im,1))*pixel_size_ax,log_im(:,:,2));colormap(gray);
xlabel('Lateral distance [mm]','FontSize',14)
ylabel('Axial distance [mm]','Fontsize',14)
axis image; title('Image 2')

im1_fft = fftshift(fft2(im(:,:,1)));
x_fft = (-0.5:(1/size(im(:,:,1),2)):(1-(1/size(im(:,:,1),1)))-0.5)';
y_fft = (-0.5:(1/size(im(:,:,1),1)):(1-(1/size(im(:,:,1),2)))-0.5)';
figure; imagesc(x_fft,y_fft,abs(im1_fft)); title('FFT of Image 1')


global champ_dep_est

%% BM PARAMETERS %%

%Estimated zone on the first image
%     ymin ymax   xmin xmax
Z_im=[ 150  2450 ;   5  190 ];

%mesh step
%     y  x
Grid=round([10 5]./[decim_y decim_x]);

%Blocks size
%  y  x
L = round([40 20]./[decim_y decim_x]);

%interpolation of search regions ([1 1] if analytical estimator)
%   y x
vI=[1 1];

%initial translation limits
%[yhaut ybas xgauche xdroite]
dP_ini_limit = [10 10 ; 0 0];

%tolerance of search regions reported to the blocks size
%[yhaut ybas xgauche xdroite]
tR=[1 1 ; 1 1];

% Local estimation type
% 1 = SAD ; 2 = normalized correlation; 3 = SAD + gradient descend search;
% 4 = elimination successive; 6 = analytical estimator;
type_correl = 6;

US_IRM = 1;



tic
%%%% estimation %%%%%%
if type_correl == 5 % corresponding to LSP estimator applied to the corr function
    [champ_dep_est x_grid y_grid]=BM(im(:,:,1),im(:,:,2),Z_im,L,Grid,vI,type_correl,dP_ini_limit,tR,f,US_IRM);
elseif type_correl == 6 % corresponding to LSP estimator applied to signals
    [champ_dep_est x_grid y_grid]=BM(im(:,:,1),im(:,:,2),Z_im,L,Grid,vI,type_correl,dP_ini_limit,tR,f,US_IRM);
else % with classical cost functions
    [champ_dep_est x_grid y_grid]=BM(im(:,:,1),im(:,:,2),Z_im,L,Grid,vI,type_correl,dP_ini_limit,tR,US_IRM);
end
estimation_time = toc

%%%% motion dense field processing %%%%%%
depl_ax = champ_dep_est(:,:,1);
depl_lat = champ_dep_est(:,:,2);
[X_grid Y_grid] = meshgrid(x_grid,y_grid);
[pts_lat pts_ax] = meshgrid(x_grid(1):x_grid(end),y_grid(1):y_grid(end));
chp_dep_int(:,:,1) = interp2(X_grid,Y_grid,depl_ax,pts_lat,pts_ax,'cubic');
chp_dep_int(:,:,2) = interp2(X_grid,Y_grid,depl_lat,pts_lat,pts_ax,'cubic');

plot_results(log_im(:,:,1),chp_dep_int(:,:,1),chp_dep_int(:,:,2),pts_ax,pts_lat,pixel_size_ax,pixel_size_lat)

% Image registration using the estimated dense motion field
% [r_corr_local,NCC_global,im1,im2_initial,im2c]=registration(im,chp_dep_int,pts_ax(1):pts_ax(end),pts_lat(1):pts_lat(end),L);
% x_corr = ((pts_lat(1)-10):L(2):(pts_lat(end)-10-2*L(2)))*pixel_size_lat;
% y_corr = ((pts_ax(1)-10):L(1):(pts_ax(end)-10-2*L(1)))*pixel_size_ax;
% figure; imagesc(x_corr,y_corr,r_corr_local*100); colorbar
% xlabel('Lateral distance [mm]','Fontsize',14)
% ylabel('Axial distance [mm]','Fontsize',14)
% title(['Moyenne = ',num2str(mean2(r_corr_local)*100,5),'///','Ecart type
% = ',num2str(std2(r_corr_local)*100,5)])
