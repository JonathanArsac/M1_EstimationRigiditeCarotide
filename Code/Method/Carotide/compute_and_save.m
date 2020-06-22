clear; close all; clc

x_min = 170;
x_max = 580;
y_min = 130;
y_max = 370;
Z_im=[ y_min  y_max ;   x_min   x_max]; % ROI
Grid=[10 10]; % distance between two pixels to estimate
L=[30 30]; % block size
vI=[1 1]; % interpolation factor of search zones
nb_it=1;
dP_ini_limit = [0 0 ; 0 0];
tR=[1 1 ;1 1]; % maximum displacement that can be estimated
type_correl = 1;
type_calcul_dep = 1;
type_interp = 2;
f = [0.13 0.13];
US_IRM = 0;

load IM_0020.mat
im = double(im);

debut = 1; pas = 1; fin = size(im,3);

out = "/home/foxiler/git/M1_EstimationRigiditeCarotide/Donnees/save_estimation/estimation_IM_0020.mat";
k=1;
progress_bar=timerbar(0,'calculating deplacement ...');
for i = debut:pas:fin-pas
    timerbar((i-debut)/(fin-debut), progress_bar);
    [chp_dep_est,x_grid,y_grid]=bdbm(im(:,:,[i, i+pas]),Z_im,L,Grid,vI,nb_it,type_correl,type_calcul_dep,dP_ini_limit,tR,f,US_IRM);
    [chp_dep_int,pts_ax,pts_lat] = denseField(chp_dep_est,x_grid, y_grid,Grid,type_interp);
    deplacement_bdbm(:,:,:,:,k) = chp_dep_est;
    deplacement_x(:,:,k) = chp_dep_int(:,:,2);
    deplacement_y(:,:,k) = chp_dep_int(:,:,1);
    k=k+1;
end
timerbar(1, progress_bar);
estimation.frames = debut:pas:fin;
estimation.chp_dep_est = deplacement_bdbm;
estimation.x_grid = x_grid;
estimation.y_grid = y_grid;
estimation.deplacement_x = deplacement_x;
estimation.deplacement_y = deplacement_y;
estimation.x_dep_grid = pts_lat;
estimation.y_dep_grid = pts_ax;
save(out, "estimation");
close(progress_bar);
