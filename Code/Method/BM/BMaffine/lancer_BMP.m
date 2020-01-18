clear; close all;
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill in the path of images location and the result file name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%schema = imread('parameters.JPG'); figure; imagesc(schema); axis image


% Name of the result file
resultat = 'result';

decim_x = 2;
decim_y = 2;

load RF2D1
im(:,:,1) = RF(1:decim_y:end,1:decim_x:end);

load RF2D2
im(:,:,2) = RF(1:decim_y:end,1:decim_x:end);

clear RF

%% IMAGE PARAMETERS %%

% normalized frequency vector
%     fy   fx
f = [0.12 0.08].*[decim_y decim_x];

% pixel size in mm
pixel_size_ax = 0.025*decim_y;
pixel_size_lat = 0.1*decim_x; %Pour l'echographe Ultrasonix taille pixel 40 mm (taille sonde) /128 lignes = 0.3125 mm

%taille de l'image
Tim=size(im(:,:,1));

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
clear im1_fft
clear log_im

global champ_dep_est

%% BM PARAMETERS %%

%Estimated zone on the first image
%     ymin ymax   xmin xmax

Z_im=[ 750  2250 ;   150  600 ];


% Z_im=[ 1500  4500 ; 400 1100 ];
%mesh step
%     y  x
Grid=round([20 10]./[decim_y decim_x]);

%Blocks size
%  y  x
L = round([91 61]./[decim_y decim_x]);
if mod(L(1),2)==0
    L(1)=L(1)+1;
end
if mod(L(2),2)==0
    L(2)=L(2)+1;
end
    %interpolation of search regions ([1 1] if analytical estimator)
%   y x
vI=[1 1];

%initial translation limits
%[yhaut ybas xgauche xdroite]
dP_ini_limit = [0 0 ; 0 0];

%tolerance of search regions reported to the blocks size
%[yhaut ybas xgauche xdroite]
% tR=[12 12 ; 4 4];
tR=[4 4 ; 1 1];
% Local estimation type
% 1 = SAD ; 2 = normalized correlation; 3 = SAD + gradient descend search;
% 4 = elimination successive; 6 = analytical estimator (translational model);9 = analytical estimator(affin model)
type_correl1 = 1; type_correl2 = 9;

if (type_correl2 == 9)||(type_correl2 == 12)
    type_fenetrage='rectwin';%choisir type d'apodistion pour l'estimation des parametre affine
else
    type_fenetrage=NaN;
end
tic
%%%% estimation %%%%%%
[champ_dep_est x_grid y_grid ]=...
        BMP(im(:,:,1),im(:,:,2),Z_im,L,Grid,vI,type_correl1,type_correl2,tR,f,type_fenetrage);
estimation_time = toc;
log_im(:,:,1) = log(abs(hilbert(im(:,:,1))));
log_im(:,:,2) = log(abs(hilbert(im(:,:,2))));

if (type_correl2 == 9)||(type_correl2 == 12)
    %fonction pour convertir les parametre affine en champ dense de
    %mouvement au niveau du pixel
    [chp_dep_int]=ConvertAffinePara(champ_dep_est,x_grid,y_grid,L,Tim,pixel_size_ax,pixel_size_lat);   
    
    chp_dep_int=chp_dep_int(y_grid(1):y_grid(end),x_grid(1):x_grid(end),:);
    
      
    [X_grid Y_grid] = meshgrid(x_grid,y_grid);
    [pts_lat pts_ax] = meshgrid(x_grid(1):x_grid(end),y_grid(1):y_grid(end));   

%  save chpdep1 chp_dep_int pts_lat pts_ax x_grid y_grid champ_dep_est
%  im_result

    plot_resultsAffine(im(:,:,1),chp_dep_int(:,:,1),chp_dep_int(:,:,2),pts_ax(:,1),pts_lat(1,:),pixel_size_ax,pixel_size_lat,champ_dep_est)

else  
    depl_ax = champ_dep_est(:,:,1);
    depl_lat = champ_dep_est(:,:,2);
    [X_grid Y_grid] = meshgrid(x_grid,y_grid);
    [pts_lat pts_ax] = meshgrid(x_grid(1):x_grid(end),y_grid(1):y_grid(end));
    chp_dep_int(:,:,1) = interp2(X_grid,Y_grid,depl_ax,pts_lat,pts_ax,'cubic');
    chp_dep_int(:,:,2) = interp2(X_grid,Y_grid,depl_lat,pts_lat,pts_ax,'cubic');

%  save chpdep1T chp_dep_int pts_lat pts_ax x_grid y_grid champ_dep_est
%  im_result2


    plot_results(im(:,:,1),chp_dep_int(:,:,1),chp_dep_int(:,:,2),pts_ax(:,1),pts_lat(1,:),pixel_size_ax,pixel_size_lat)

end

 



 
 
 
 