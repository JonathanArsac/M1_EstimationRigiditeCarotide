clear; close all; clc

Z_im=[ 250  3000 ;   7  293 ];
Grid=[15 4];
L=[40 20];
vI=[1 1];
nb_it = 1;
tR=[15 0; 1 1];
dP_ini_limit = [0 0 ; 0 0];
type_correl = 6;
f = [0.1 0.06];
type_calcul_dep= 1;
type_interp=1;
US_IRM = 0;

load imUScardio
im = seq; clear seq;

debut = 1; pas = 1; fin = 30;

for i = debut:pas:fin
   
    image(:,:,1)=im(:,:,i); image(:,:,2)=im(:,:,i+pas);
    
    [chp_dep_est x_grid y_grid]=bdbm2(image,Z_im,L,Grid,vI,nb_it,type_correl,type_calcul_dep,dP_ini_limit,tR,f,US_IRM);   
    [chp_dep_int,pts_ax,pts_lat] = denseField(chp_dep_est,x_grid, y_grid,Grid,type_interp);

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




