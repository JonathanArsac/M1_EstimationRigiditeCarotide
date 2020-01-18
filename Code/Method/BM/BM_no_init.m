function [champ_dep_est, x_grid, y_grid]=BM_no_init(im1,im2,Z_im,L,Grid,vI,type_correl,dP_ini_limit,tR,f,US_IRM);

%2-D motion estimation with block matching (speckle tracking)

%[chp_dep_est, x_grid, y_grid]=speckle_tracking_simple(im1,im2,Z_im,L,[Grid],[vI],[t_cor],[dP_i],[tR]);

%INPUT
%  im1 :   1st image
%  im2 :   2nd image
%  Z_im :  estimated zone [ymin ymax xmin xmax]
%  L :     block size [Lv Lu]
%  Grid :  mesh grid on the 1st image with points to be estimated [Gy Gx]
%  vI :    interpolation factors of the 2nd image [int_y int_x]
%  type_cor : type of local 2D spatial delay estimation.
%                                           1 = SAD ; 2 = normalized correlation; 3 = SAD + gradient descend search;
%                                           4 = elimination successive; 5 = LSP on correlation function;
%                                           6 = LSP on signals
%
%  dP_ini_limit :  limits of the initial displacement [yup ydown xleft xright]
%  tR :    tolerance of search regions reported to the blocks size [yup ydown xleft xright]
%  f = [fy fx] : frequencies along y and x directions
%  US_IRM : 1 if US, 0 if IRM
%OUTPUT
%  chp_dep_est : estimated motion 3D
%     dim 1 : ligne (y)
%     dim 2 : colonne (x)
%     dim 3 : direction (1=y, 2=x)
%  x_grid : lateral grid of estimated pixels
%  y_grid : axial grid of estimated pixels

% Adrian Basarab, Philippe Delachartre, 2006



if ~exist('Grid','var')
    Grid = L;
end
if ~exist('vI','var')
    vI = [3 3];
end
if ~exist('type_correl','var')
    type_correl = 1;
end
if ~exist('dP_ini','var')
    dP_ini = [0 0];
end
if ~exist('tR','var')
    tR=[1  1 ; 1 1];;
end
if ~exist('f','var')
    f=[0 0];
end

if type_correl == 5
    vI = [1 1]; % no interpolation in the case of LSE phase estimator
    tR = [1 1 ; 1 1];
elseif type_correl == 6
    vI = [1 1]; % no interpolation in the case of LSE phase estimator
    tR = [1 1 ; 1 1];
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
centre_im_N=fix(N_max/2)*Grid(2)+debut;     %sur l'image
N_mat_dep=fix(N_max/2)+1;                   %dans la matrice

%pixel central de départ
P_dep=[fix(L(1)/2)+tR(1,1)+Z_im(1,1) centre_im_N];

%grille de calcul
x_grid=debut:Grid(2):debut+(N_max-1)*Grid(2);
y_grid=P_dep(1):Grid(1):P_dep(1)+M_max*Grid(1);
%nombre de points d'estimation
nb_pts_est=(N_max-1)*M_max;

%initialisation nombres de points calculés
nb_pts_calcul=0;

%nombres de boucles de calcul
nb_boucle=N_max+M_max-N_mat_dep+1;

%initialisation matrice champ de déplacement
champ_dep_est=NaN*ones(M_max,N_max,2);

%grille d'interpolation
[ix,iy]=meshgrid(1:tZ(2)+1,1:tZ(1)+1);
[ixi,iyi]=meshgrid(1:1/vI(2):tZ(2)+1,1:1/vI(1):tZ(1)+1);


% % estimate the initial displacement of the first block
% B1 = im1(P_dep(1)-ceil(L(1)/2)+1:P_dep(1)-ceil(L(1)/2)+L(1),...
%     P_dep(2)-ceil(L(2)/2)+1:P_dep(2)-ceil(L(2)/2)+L(2));
% ZRtest = im2(P_dep(1)-ceil(L(1)/2)+1-dP_ini_limit(1,1):P_dep(1)-ceil(L(1)/2)+L(1)+dP_ini_limit(1,2),...
%     P_dep(2)-ceil(L(2)/2)+1-dP_ini_limit(2,1):P_dep(2)-ceil(L(2)/2)+L(2)+dP_ini_limit(2,2));
% 
% [pos_best_corr, vcorr]=CostFunctions(B1,ZRtest,[1 1],1);

dP_ini=[0 0];%(pos_best_corr-1-dP_ini_limit(:,1)');


%progress bar
progress_bar=timerbar(0,'calculating...');

%parcours de l'image en triangle
for cercle_it=0:nb_boucle
    %bornes de calcul en largeur
    borne_min_N=max(1-N_mat_dep,-cercle_it);
    borne_max_N=min(N_max-N_mat_dep,cercle_it);
    %progress bar - estimation du temps restant
    timerbar(nb_pts_calcul/nb_pts_est,progress_bar);
    nb_pts_calcul=nb_pts_calcul+(borne_max_N-borne_min_N+1);
    %points du triangle
    for jN=borne_min_N:borne_max_N
        jM=min(cercle_it-abs(jN),M_max);
        clc
        %déplacement supposé (répercution du déplacement trouvé sur les pixels voisins)
%         if jM==0 & jN==0
%             dP=dP_ini;
%         else
%             if jN~=0
%                 dP(1)=round(champ_dep_est(jM+1,N_mat_dep+jN-sign(jN),1));
%             else
%                 dP(1)=round(champ_dep_est(jM,N_mat_dep,1));
%             end
%             if jM~=0
%                 dP(2)=round(champ_dep_est(jM,N_mat_dep+jN,2));
%             else
%                 dP(2)=round(champ_dep_est(1,N_mat_dep+jN-sign(jN),2));
%             end
%         end
        dP = [0 0];
        %point central de la zone d'étude
        P=P_dep+Grid.*[jM jN];
        %pattern initial
        Pat_ini=im1(P(1)-ceil(L(1)/2)+1:P(1)-ceil(L(1)/2)+L(1),P(2)-ceil(L(2)/2)+1:P(2)-ceil(L(2)/2)+L(2));
        
        if type_correl==5
            ZR = im2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2));
            
            [dPx dPy] = PhaseAdjust(Pat_ini,ZR,f(2),f(1));
            dep = [dPy dPx] + dP;
            
        elseif type_correl==6
            if US_IRM
                ZR = im2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2));
                
                [dPx dPy] = PhaseAdjust3(P1_im1(P(1)-ceil(L(1)/2)+1:P(1)-ceil(L(1)/2)+L(1),P(2)-ceil(L(2)/2)+1:P(2)-ceil(L(2)/2)+L(2)),...
                    P2_im1(P(1)-ceil(L(1)/2)+1:P(1)-ceil(L(1)/2)+L(1),P(2)-ceil(L(2)/2)+1:P(2)-ceil(L(2)/2)+L(2))...
                    ,P1_im2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2))...
                    ,P2_im2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2)),f(2),f(1));
                %             dPx = 0;
                dep = [dPy dPx] + dP;
            else
                
                [dPx dPy] = PhaseAdjust2(im1(P(1)-ceil(L(1)/2)+1:P(1)-ceil(L(1)/2)+L(1),P(2)-ceil(L(2)/2)+1:P(2)-ceil(L(2)/2)+L(2)),...
                    im2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2)),...
                    f(2),f(1),US_IRM);
                dep = [dPy dPx] + dP;
            end
            
        else
            %zone de recherche
            PsZ=P-ceil(L/2)+1-tR(:,1)'+dP;     %coordonnée pixel haut-gauche de la zone de recherche
            ZR=im2(PsZ(1):PsZ(1)+tZ(1),PsZ(2):PsZ(2)+tZ(2));  %zone de recherche
            %interpolation
            Pat_ZR=interp2(ix,iy,ZR,ixi,iyi);
            %estimation du champ de déplacement en Pp
            
            [pos_best_corr, vcorr]=CostFunctions(Pat_ini,Pat_ZR,vI,type_correl);
            
            %déplacement réel (ajout du déplacement de départ)
            dep=((pos_best_corr-1)./vI-tR(:,1)')+dP;
        end
        %construction de la matrice du champ de déplacement
        champ_dep_est(jM+1,N_mat_dep+jN,:)=dep;
    end
end
close(progress_bar);

