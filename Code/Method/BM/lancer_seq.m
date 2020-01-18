clear; close all; clc

Z_im=[ 100  1850 ;   20  100 ];
Grid=[2 1];
L=[10 5];
vI=[1 1];
dP_ini_limit = [5 5 ; 0 0];
tR=[1 1 ;1 1];
type_correl = 6;
f = [0.12 0.1];

debut = 5; pas = 1; fin = 20;
x=5:70;
y=120:1700;

load imUS

for i = debut:pas:fin
   
    [champ_dep_est x_grid y_grid]=fastBM(im(:,:,i),im(:,:,i+pas),Z_im,L,Grid,vI,type_correl,dP_ini_limit,tR,f);
    
    depl_ax = champ_dep_est(:,:,1);
    depl_lat = champ_dep_est(:,:,2);
    [X_grid Y_grid] = meshgrid(x_grid,y_grid);
    [pts_lat pts_ax] = meshgrid(x_grid(1):x_grid(end),y_grid(1):y_grid(end));
    chp_dep_int(:,:,1) = interp2(X_grid,Y_grid,depl_ax,pts_lat,pts_ax,'cubic');
    chp_dep_int(:,:,2) = interp2(X_grid,Y_grid,depl_lat,pts_lat,pts_ax,'cubic');

    if i==debut
        depl_lat_cumule(:,:,1) = chp_dep_int(:,:,2);
        depl_ax_cumule(:,:,1) = chp_dep_int(:,:,1);
    else
        [depl_lat_cumule,depl_ax_cumule] = ...
            displacementReg(pts_ax,pts_lat,chp_dep_int(:,:,1),chp_dep_int(:,:,2),depl_lat_cumule,depl_ax_cumule);
    end

    axial_motion = sum(depl_ax_cumule(y,x,1:(i-debut+pas)/pas),3);
    Nc = floor(0.1*[size(axial_motion,1)])+1;
    eyy = derivate2(axial_motion,Nc)*100;
    imagesc(x,((y(1)+floor(Nc/2)):((y(end)-floor(Nc/2)))),-eyy); colorbar; %colormap(gray)
    H((i-debut+pas)/pas) = getframe
end



