clear; close all; clc

x_min = 150;
x_max = 550;
y_min = 100;
y_max = 450;
Z_im=[ y_min  y_max ;   x_min   x_max]; % ROI
Grid=[10 10]; % distance between two pixels to estimate
L=[30 30]; % block size
vI=[1 1]; % interpolation factor of search zones
nb_it=1;
dP_ini_limit = [0 0 ; 0 0];
tR=[1 1 ;1 1]; % maximum displacement that can be estimated
type_correl = 1;
type_calcul_dep = 1;
type_interp = 1;
f = [0.13 0.13];
US_IRM = 0;

load IM_0021.mat
im = double(im);

debut = 1; pas = 10; fin = 100;

% For motion vector plot
pas_x = 10; pas_y = 10;

out = "/home/foxiler/git/M1_EstimationRigiditeCarotide/Donnees/save_estimation/IM_0021.mat";
k=1;
progress_bar=timerbar(0,'calculating ...');
for i = debut:pas:fin
    %% 
    timerbar(0);
    timerbar((i-debut) / (fin-debut));
    image(:,:,1)=im(:,:,i); image(:,:,2)=im(:,:,i+pas);
    
    [chp_dep_est,x_grid,y_grid]=bdbm(image,Z_im,L,Grid,vI,nb_it,type_correl,type_calcul_dep,dP_ini_limit,tR,f,US_IRM);   
    [chp_dep_int,pts_ax,pts_lat] = denseField(chp_dep_est,x_grid, y_grid,Grid,type_interp);
    estimation.chp_dep_est(k,:,:,:,:,:) = chp_dep_est;
    if i==debut
        depl_lat_cumule(:,:,1) = chp_dep_int(:,:,2);
        depl_ax_cumule(:,:,1) = chp_dep_int(:,:,1);
        estimation.x_grid = x_grid;
        estimation.y_grid = y_grid;
        estimation.pts_ax = pts_ax;
        estimation.pts_lat = pts_lat;
    else
        [depl_lat_cumule,depl_ax_cumule] = ...
            displacementReg(pts_ax,pts_lat,chp_dep_int(:,:,1),chp_dep_int(:,:,2),depl_lat_cumule,depl_ax_cumule);
    end
    k=k+1;
    
end
close(progress_bar);
estimation.depl_lat_cumule = depl_lat_cumule;
estimation.depl_ax_cumule = depl_ax_cumule;
save(out, "estimation");
