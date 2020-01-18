function [champ_dep_est, x_grid, y_grid]=BM2(im1,im2,Z_im,L,Grid,vI,type_cor1,type_cor2,tR,f);

% 2-D motion estimation with block matching (speckle tracking)
% Designed to estimate large local motion
%   - classical scan of nodes (from top to bottom and from left to right)
%   - no initialisation using the results of neighbouring nodes
%   - for each node:
%           - integer pixel accuracy search using a tolerance of tR
%           - subpixel search around the position found at previous step:
%                   - using LSP estimator OR
%                   - using classical cost functions

%[chp_dep_est, x_grid, y_grid]=speckle_tracking_simple(im1,im2,Z_im,L,[Grid],[vI],[t_cor],[dP_i],[tR]);

%INPUT
%  im1 :   1st image
%  im2 :   2nd image
%  Z_im :  estimated zone [ymin ymax xmin xmax]
%  L :     block size [Lv Lu]
%  Grid :  mesh grid on the 1st image with points to be estimated [Gy Gx]
%  vI :    interpolation factors of the 2nd image [int_y int_x]
%  type_cor1 and type_cor2 : type of local 2D spatial delay estimation.
%                                           1 = SAD ; 2 = normalized correlation; 3 = SAD + gradient descend search;
%                                           4 = elimination successive; 5 = LSP on correlation function;
%                                           6 = LSP on signals
%
% type_cor1 for the integer pixel research
%
%  tR :    tolerance of search regions reported to the blocks size [yup ydown xleft xright]
%  f = [fy fx] : frequencies along y and x directions
%OUTPUT
%  chp_dep_est : estimated motion 3D
%     dim 1 : ligne (y)
%     dim 2 : colonne (x)
%     dim 3 : direction (1=y, 2=x)
%  x_grid : lateral grid of estimated pixels
%  y_grid : axial grid of estimated pixels

% Adrian Basarab, Philippe Delachartre, 2009



if ~exist('Grid','var')
    Grid = L;
end
if ~exist('vI','var')
    vI = [3 3];
end
if ~exist('type_cor2','var')
    type_cor2 = 1;
end
if ~exist('tR','var')
    tR=[1  1 ; 1 1];
end
if ~exist('f','var')
    f=[0 0];
end

if type_cor2 == 5
    vI = [1 1]; % no interpolation in the case of LSE phase estimator
    tR = [1 1 ; 1 1];
elseif type_cor2 == 6
    vI = [1 1]; % no interpolation in the case of LSE phase estimator
    % compute the analytical signals for im1 and im2

    [R1,R2]=ashahn(im1);
    [S1,S2]=ashahn(im2);
    P1_im1=atan2(imag(R1),real(R1));
    P2_im1=atan2(imag(R2),real(R2));
    P1_im2=atan2(imag(S1),real(S1));
    P2_im2=atan2(imag(S2),real(S2));

else
    f = [0 0]; % no frequency vector needed in the case of classical estimators
end

%taille d'un pattern de recherche de la zone de recherche
tZ=L-1+sum(tR');

%taille de la zone d'étude
long_Z_im=diff(Z_im(1,:))-tZ(1);
M_max=fix(long_Z_im/Grid(1))+1;
larg_Z_im=diff(Z_im(2,:))-tZ(2);
N_max=fix(larg_Z_im/Grid(2));
%point de calcul le plus à gauche
debut=tZ(2)-tR(2,2)+Z_im(2,1)+fix((larg_Z_im-(N_max-1)*Grid(2))/2);
%centre de départ
left_im_N=debut;             %sur l'image
N_mat_dep=1;                 %dans la matrice

%pixel central de départ (noeud en haut à gauche)
P_dep=[fix(L(1)/2)+tR(1,1)+Z_im(1,1) left_im_N];

%grille de calcul
x_grid=debut:Grid(2):debut+(N_max-1)*Grid(2);
y_grid=P_dep(1):Grid(1):P_dep(1)+M_max*Grid(1);

%initialisation matrice champ de déplacement
champ_dep_est=NaN*ones(M_max+1,N_max,2);

%grille d'interpolation
[ix,iy]=meshgrid(1:tZ(2)+1,1:tZ(1)+1);
[ixi,iyi]=meshgrid(1:1/vI(2):tZ(2)+1,1:1/vI(1):tZ(1)+1);

%progress bar
progress_bar=timerbar(0,'calculating...');
nb_pts_est=N_max*(M_max+1);
%initialisation nombres de points calculés
nb_pts_calcul=0;

% nodes scan
for jM=0:M_max

    nb_pts_calcul = jM+nb_pts_calcul;
    timerbar(nb_pts_calcul/nb_pts_est,progress_bar);

    for jN=0:(N_max-1)
        clc

        % Integer pixel search
        P=P_dep+Grid.*[jM jN];
        Pat_ini=im1(P(1)-ceil(L(1)/2)+1:P(1)-ceil(L(1)/2)+L(1),P(2)-ceil(L(2)/2)+1:P(2)-ceil(L(2)/2)+L(2));
        PsZ=P-ceil(L/2)+1-tR(:,1)';     %coordonnée pixel haut-gauche de la zone de recherche
        ZR=im2(PsZ(1):PsZ(1)+tZ(1),PsZ(2):PsZ(2)+tZ(2));  %zone de recherche
        %estimation du champ de déplacement en Pp
        [pos_best_corr, vcorr]=CostFunctions(Pat_ini,ZR,[1 1],type_cor1);
        %déplacement réel (ajout du déplacement de départ)
        dep_init=((pos_best_corr-1)./[1 1]-tR(:,1)');

        if max(max(abs(Pat_ini)))==0
            dP = [0 0];
        else
            dP=dep_init;
        end

        if max(max(abs(Pat_ini)))==0
            dep = [0 0];
        else
            if type_cor2==5
                ZR = im2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2));
                [dPx dPy] = PhaseAdjust(Pat_ini,ZR,f(2),f(1));
                dep = [dPy dPx] + dP;

            elseif type_cor2==6

                [dPx dPy] = PhaseAdjust3(P1_im1(P(1)-ceil(L(1)/2)+1:P(1)-ceil(L(1)/2)+L(1),P(2)-ceil(L(2)/2)+1:P(2)-ceil(L(2)/2)+L(2)),...
                    P2_im1(P(1)-ceil(L(1)/2)+1:P(1)-ceil(L(1)/2)+L(1),P(2)-ceil(L(2)/2)+1:P(2)-ceil(L(2)/2)+L(2))...
                    ,P1_im2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2))...
                    ,P2_im2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2)),f(2),f(1));
                dep = [dPy dPx] + dP;

            else
                %zone de recherche
                PsZ=P-ceil(L/2)+1-[1 1]+dP;     %coordonnée pixel haut-gauche de la zone de recherche
                ZR=im2(PsZ(1):PsZ(1)+tZ(1),PsZ(2):PsZ(2)+tZ(2));  %zone de recherche
                %interpolation
                Pat_ZR=interp2(ix,iy,ZR,ixi,iyi);
                %estimation du champ de déplacement en Pp

                [pos_best_corr, vcorr]=CostFunctions(Pat_ini,Pat_ZR,vI,type_cor2);

                %déplacement réel (ajout du déplacement de départ)
                dep=((pos_best_corr-1)./vI-[1 1])+dP;
            end
        end
        %construction de la matrice du champ de déplacement
        champ_dep_est(jM+1,jN+1,:)=dep;
    end
end
close(progress_bar);

