clear; close all; clc

x_min = 170;
x_max = 580;
y_min = 130;
y_max = 370;
Z_im=[ y_min  y_max ;   x_min   x_max]; % ROI
Grid=[10 10]; % distance between two pixels to estimate
L=[30 30]; % block size
vI=[1 1]; % interpolation factor of search zones
nb_it=2;
dP_ini_limit = [0 0 ; 0 0];
tR=[1 1 ;1 1]; % maximum displacement that can be estimated
type_correl = 1;
type_calcul_dep = 1;
type_interp = 1;
f = [0.13 0.13];
US_IRM = 0;

load IM_0020.mat
im = double(im);

% debut = 1; pas = 1; fin = size(im,3);
debut = 1; pas = 1; fin = size(im,3);

out = "/home/foxiler/git/M1_EstimationRigiditeCarotide/Donnees/save_estimation/estimation_IM_0020.mat";
[chp_dep_est,x_grid,y_grid]=bdbm(im(:,:,debut:pas:fin),Z_im,L,Grid,vI,nb_it,type_correl,type_calcul_dep,dP_ini_limit,tR,f,US_IRM);
k=1;
progress_bar=timerbar(0,'calculating deplacement ...');
for i = debut:pas:fin-1
    timerbar(k/(fin-debut), progress_bar);
    [chp_dep_int,pts_ax,pts_lat] = denseField(chp_dep_est(:,:,:,:,k),x_grid, y_grid,Grid,type_interp);
    deplacement_x(:,:,k) = chp_dep_int(:,:,2);
    deplacement_y(:,:,k) = chp_dep_int(:,:,1);
    k=k+1;
end
close(progress_bar);
estimation.x_dep_grid = pts_lat;
estimation.y_dep_grid = pts_ax;
estimation.deplacement_x = deplacement_x;
estimation.deplacement_y = deplacement_y;
estimation.chp_dep_est = chp_dep_est;
estimation.x_grid = x_grid;
estimation.y_grid = y_grid;
estimation.frames = debut:pas:fin;
save(out, "estimation");
