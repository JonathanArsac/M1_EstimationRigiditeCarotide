clear; close all; clc

Z_im=[ 50  160 ;   10   140]; % ROI
Grid=[3 3]; % distance between two pixels to estimate
L=[10 10]; % block size
vI=[1 1]; % interpolation factor of search zones
nb_it=1;
dP_ini_limit = [0 0 ; 0 0];
tR=[1 1 ;1 1]; % maximum displacement that can be estimated
type_correl = 6;
type_calcul_dep = 1;
type_interp = 1;
f = [0.13 0.13];
US_IRM = 0;

load imTaggedMRI.mat
im = double(im(50:250,1:160,:));

debut = 1; pas = 1; fin = 11;

% For motion vector plot
pas_x = 5; pas_y = 5;


for i = debut:pas:fin
   
    image(:,:,1)=im(:,:,i); image(:,:,2)=im(:,:,i+pas);
    
    [chp_dep_est x_grid y_grid]=bdbm(image,Z_im,L,Grid,vI,nb_it,type_correl,type_calcul_dep,dP_ini_limit,tR,f,US_IRM);   
    [chp_dep_int,pts_ax,pts_lat] = denseField(chp_dep_est,x_grid, y_grid,Grid,type_interp);

    if i==debut
        depl_lat_cumule(:,:,1) = chp_dep_int(:,:,2);
        depl_ax_cumule(:,:,1) = chp_dep_int(:,:,1);
    else
        [depl_lat_cumule,depl_ax_cumule] = ...
            displacementReg(pts_ax,pts_lat,chp_dep_int(:,:,1),chp_dep_int(:,:,2),depl_lat_cumule,depl_ax_cumule);
    end

    axial_motion(:,:,(i-debut+pas)/pas) = sum(depl_ax_cumule(:,:,1:(i-debut+pas)/pas),3);
    lateral_motion(:,:,(i-debut+pas)/pas) = sum(depl_lat_cumule(:,:,1:(i-debut+pas)/pas),3);
    
    % Plot motion vectors
    x=pts_lat(1):pts_lat(end);
    y=pts_ax(1):pts_ax(end);
    depl_ax_quiver = chp_dep_int(1:pas_y:end,1:pas_x:end,1);
    depl_lat_quiver = chp_dep_int(1:pas_y:end,1:pas_x:end,2);
    [ax_x_quiver ax_y_quiver] = meshgrid(floor(x),floor(y));
    image_irm = image(floor(y),floor(x),1);
    imagesc(floor(x),floor(y),image_irm); colormap(gray); %axis image;
    hold on
    quiver(floor(ax_x_quiver(1:pas_y:end,1:pas_x:end)),floor(ax_y_quiver(1:pas_y:end,1:pas_x:end))...
        ,5*depl_lat_quiver,5*depl_ax_quiver,0,'Color','w');colormap(gray);
    H(i) = getframe
    
end

figure;
J1 = DeformMesh(im(40:160,20:120,:),[5 5],18:90,28:100,axial_motion,lateral_motion,x_grid-20,y_grid-40,[2 2])
movie(J1,100,5)
% J1 = DeformMesh45d(im(100:200,20:110,:),[3 3],10:90,20:90,axial_motion,lateral_motion,x_grid-20,y_grid-100,[2 2])



