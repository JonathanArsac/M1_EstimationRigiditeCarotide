function [champ_dep_est, x_grid, y_grid]=BM_monogenic(im1,im2,Z_im,L,Grid,dP_ini_limit,tR);

%2-D motion estimation with block matching using monogenic signal


%INPUT
%  im1 :   1st image
%  im2 :   2nd image
%  Z_im :  estimated zone [ymin ymax xmin xmax] 
%  L :     block size [Lv Lu]
%  Grid :  mesh grid on the 1st image with points to be estimated [Gy Gx]
%                                           
%  dP_ini_limit :  limits of the initial displacement [yup ydown xleft xright]
%  tR :    tolerance of search regions reported to the blocks size [yup
%  ydown xleft xright]
%OUTPUT
%  chp_dep_est : estimated motion 3D
%     dim 1 : ligne (y)
%     dim 2 : colonne (x)
%     dim 3 : direction (1=y, 2=x)
%  x_grid : lateral grid of estimated pixels
%  y_grid : axial grid of estimated pixels

% Adrian Basarab, Philippe Delachartre, 2006
% Tanguy Maltaverne, 2009



[a1,phi1,theta1,freq1] = monogenic2(im1);
[a2,phi2,theta2,freq2] = monogenic2(im2);

% [f, h1f, h2f, amp1, t1, p1, freq1] = ...
%     monofilt(im1, 3, 4, 2, 0.65, 0);
% [f, h1f, h2f, amp2, t2, p2, freq2] = ...
%     monofilt(im2, 3, 4, 2, 0.65, 0);
% a1 = amp1{3}; a2 = amp2{3};
% theta1 = t1{3}; theta2 = t2{3};
% phi1 = p1{3}; phi2 = p2{3};

% figure; imagesc(theta1); title('Local orientation for image 1');
% figure; imagesc(theta2); title('Local orientation for image 2');
% figure; imagesc(phi1); title('Local phase for image 1');
% figure; imagesc(phi2); title('Local phase for image 2');

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


% estimate the initial displacement of the first block
B1 = im1(P_dep(1)-ceil(L(1)/2)+1:P_dep(1)-ceil(L(1)/2)+L(1),...
    P_dep(2)-ceil(L(2)/2)+1:P_dep(2)-ceil(L(2)/2)+L(2));
ZRtest = im2(P_dep(1)-ceil(L(1)/2)+1-dP_ini_limit(1,1):P_dep(1)-ceil(L(1)/2)+L(1)+dP_ini_limit(1,2),...
    P_dep(2)-ceil(L(2)/2)+1-dP_ini_limit(2,1):P_dep(2)-ceil(L(2)/2)+L(2)+dP_ini_limit(2,2));

[pos_best_corr, vcorr]=CostFunctions(B1,ZRtest,[1 1],1);

dP_ini=(pos_best_corr-1-dP_ini_limit(:,1)');


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
%         clc
        %déplacement supposé (répercution du déplacement trouvé sur les pixels voisins)
        if jM==0 & jN==0
            dP=dP_ini;
        else
            if jN~=0
                dP(1)=round(champ_dep_est(jM+1,N_mat_dep+jN-sign(jN),1));
            else
                dP(1)=round(champ_dep_est(jM,N_mat_dep,1));
            end
            if jM~=0
                dP(2)=round(champ_dep_est(jM,N_mat_dep+jN,2));
            else
                dP(2)=round(champ_dep_est(1,N_mat_dep+jN-sign(jN),2));
            end
        end
        
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        dP = [0 0];
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        %point central de la zone d'étude
        P=P_dep+Grid.*[jM jN];  
        %pattern initial
        Pat_ini=im1(P(1)-ceil(L(1)/2)+1:P(1)-ceil(L(1)/2)+L(1),P(2)-ceil(L(2)/2)+1:P(2)-ceil(L(2)/2)+L(2));
            
        ZR = im2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2));

        [dPx dPy] = MonogenicAdjust(phi1(P(1)-ceil(L(1)/2)+1:P(1)-ceil(L(1)/2)+L(1),P(2)-ceil(L(2)/2)+1:P(2)-ceil(L(2)/2)+L(2)),...
            theta1(P(1)-ceil(L(1)/2)+1:P(1)-ceil(L(1)/2)+L(1),P(2)-ceil(L(2)/2)+1:P(2)-ceil(L(2)/2)+L(2))...
            ,freq1(P(1)-ceil(L(1)/2)+1:P(1)-ceil(L(1)/2)+L(1),P(2)-ceil(L(2)/2)+1:P(2)-ceil(L(2)/2)+L(2))...
            ,phi2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2))...
            ,theta2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2))...
            ,freq2(P(1)-ceil(L(1)/2)+1+dP(1):P(1)-ceil(L(1)/2)+L(1)+dP(1),P(2)-ceil(L(2)/2)+1+dP(2):P(2)-ceil(L(2)/2)+L(2)+dP(2)));

        dep = [dPy dPx] + dP;

        %construction de la matrice du champ de déplacement
        champ_dep_est(jM+1,N_mat_dep+jN,:)=dep;
    end
end

close(progress_bar);

