clear; close all; clc

%Z_im=[ 250  3000 ;   7  293 ];
Z_im=[ 230  3000 ;   2  300 ];
Grid=[15 4];
L=[40 20];
vI=[1 1];
tR=[15 0; 1 1];
type_correl1 = 1; type_correl2 = 6;
f = [0.1 0.06];

debut = 1; pas = 1; fin = 30;

load imUScardio
im = seq; clear seq;

for i = debut:pas:fin
   
    [champ_dep_est x_grid y_grid]=...
        BM2(im(:,:,1),im(:,:,2),Z_im,L,Grid,vI,type_correl1,type_correl2,tR,f);
   
    depl_ax = champ_dep_est(:,:,1);
    depl_lat = champ_dep_est(:,:,2);
    [X_grid Y_grid] = meshgrid(x_grid,y_grid);
    [pts_lat pts_ax] = meshgrid(x_grid(1):x_grid(end),y_grid(1):y_grid(end));
    chp_dep_int(:,:,1) = interp2(X_grid,Y_grid,depl_ax,pts_lat,pts_ax,'cubic');
    chp_dep_int(:,:,2) = interp2(X_grid,Y_grid,depl_lat,pts_lat,pts_ax,'cubic');

    if i==debut
        depl_lat_cumule(:,:,1) = chp_dep_int(:,:,2);
        depl_ax_cumule(:,:,1) = chp_dep_int(:,:,1);
        axial_motion = sum(depl_ax_cumule(:,:,1:(i-debut+pas)/pas),3);
        lateral_motion = sum(depl_lat_cumule(:,:,1:(i-debut+pas)/pas),3);
    else
        [depl_lat_cumule,depl_ax_cumule] = ...
            displacementReg(pts_ax,pts_lat,chp_dep_int(:,:,1),chp_dep_int(:,:,2),depl_lat_cumule,depl_ax_cumule);
        axial_motion = sum(depl_ax_cumule(:,:,1:(i-debut+pas)/pas),3);
        lateral_motion = sum(depl_lat_cumule(:,:,1:(i-debut+pas)/pas),3);
        depl_lat_cumule(:,:,1) = lateral_motion;
        depl_ax_cumule(:,:,1) = axial_motion;
    end

%     axial_motion = sum(depl_ax_cumule(:,:,1:(i-debut+pas)/pas),3);
    imagesc(axial_motion)
    G((i-debut+pas)/pas) = getframe
%     lateral_motion = sum(depl_lat_cumule(:,:,1:(i-debut+pas)/pas),3);
    imagesc(lateral_motion)
    H((i-debut+pas)/pas) = getframe
end



