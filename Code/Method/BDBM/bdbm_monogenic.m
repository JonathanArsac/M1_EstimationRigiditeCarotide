function [champ_dep_est, x_grid, y_grid]=bdbm_monogenic(im,Z_im,L,Grid,vI,nb_it,t_cor,type_calcul_dep,dP_ini_limit,tR,f,US_IRM);

% 2D parametric block matching motion estimation based on bilinear motion
% model

%INPUT
%  im :    images im(:,:,:)
%  Z_im :  estimated zone [ymin ymax xmin xmax]
%  L :     ROIs size [Lv Lu]
%  Grid :  mesh grid on the initial image [Gy Gx]
%  vI :    interpolation factor for each iteration [int_y int_x]
%  nb_it : number of iterations
%  t_cor : type of local 2D spatial delay estimation.
%                                           1 = SAD ; 2 = normalized correlation; 3 = SAD + gradient descend search;
%                                           6 = LSP on signals
%
%  type_calcul_dep: 1= mean,2 = median
%  dP_ini_limit :  limits of the initial displacement [yup ydown xleft xright]
%  tR :   tolerance of search regions reported to the blocks size [yhaut ybas xgauche xdroite]
%  f=[fy fx] : 2-D frequency vector
%  US_IRM : 1 if US, 0 if IRM

%OUTPUT
%  chp_dep_est : estimated bilinear parameters 5D
%     dim 1 : line (y)
%     dim 2 : colonne (x)
%     dim 3 : parameter of bilinear model (1=a, 2=b, 3=c, 4=d)
%     dim 4 : direction (1=y, 2=x)
%     dim 5 : k -> estimation between images k and k+1
%  x_grid : lateral grid of estimated pixels
%  y_grid : axial grid of estimated pixels

% Adrian Basarab, Philippe Delachartre, 2006

% Reference
%   [1] A. Basarab, H. Liebgott, F. Morestin, A. Lyshchik, T. Higashi,
%       R. Asato, P. Delachartre, A method for vector displacement estimation with
%       ultrasound imaging and its application for thyroid nodular disease,
%       Elsevier, Medical Image Analysis, 2007.


[a1,phi1,theta1,freq1] = monogenic2(im(:,:,1));
[a2,phi2,theta2,freq2] = monogenic2(im(:,:,2));


a=Grid;
if ~exist('a','var')
    Grid = L;
end
a=vI;
if ~exist('a','var')
    vI = [3 3];
end
a=nb_it;
if ~exist('a','var')
    nb_it = 2;
end
a=t_cor;
if ~exist('a','var')
    t_cor = 1;
end
a=dP_ini_limit;
if ~exist('a','var')
    dP_ini_limit = [0 0;0 0;0 0;0 0];
end
a=tR;
if ~exist('a','var')
    tR=[1  1 ; 1 1];
end

if t_cor == 6
    vI = [1 1]; % no interpolation in the case of LSE phase estimator
    tR = [1 1 ; 1 1];
    
    [R1,R2]=ashahn(im(:,:,1));
    [S1,S2]=ashahn(im(:,:,2));
    P1_im1=atan2(imag(R1),real(R1));
    P2_im1=atan2(imag(R2),real(R2));
    P1_im2=atan2(imag(S1),real(S1));
    P2_im2=atan2(imag(S2),real(S2));
    
else
    f = [0 0]; % no frequency vector needed in the case of classical estimators
end


%Numbers of images in the sequence
N = size(im,3);

global champ_dep_est
global x_grid
global y_grid

%Search zone size
tZ=L-1+sum(tR');
%Study zone size
tROI=2*L-2+sum(tR');

%taille de la zone d'étude
long_Z_im=diff(Z_im(1,:))-tROI(1);
M_max=fix(long_Z_im/Grid(1))+1;
larg_Z_im=diff(Z_im(2,:))-tROI(2);
N_max=fix(larg_Z_im/Grid(2));
%point de calcul le plus à gauche
debut=tZ(2)-tR(2,2)+Z_im(2,1)+fix((larg_Z_im-(N_max-1)*Grid(2))/2);
%centre de départ
centre_im_N=fix(N_max/2)*Grid(2)+debut;     %image coordinates
N_mat_dep=fix(N_max/2)+1;                   %coordinates in the matrix

%first node to be estimated
P_dep=[L(1)+tR(1,1)+Z_im(1,1) centre_im_N];

%grid of the nodes to be estimated
x_grid=debut:Grid(2):debut+(N_max-1)*Grid(2);
y_grid=P_dep(1):Grid(1):P_dep(1)+M_max*Grid(1);

%number of nodes to be estimated
nb_pts_est=(N_max-1)*M_max;


nb_pts_calcul=0;

%number of loops
nb_boucle=N_max+M_max-N_mat_dep+1;

champ_dep_est=NaN*ones(M_max,N_max,4,2);

% estimate the initial displacement of the first block
B1 = im(P_dep(1)-ceil(L(1)/2)+1:P_dep(1)-ceil(L(1)/2)+L(1),...
    P_dep(2)-ceil(L(2)/2)+1:P_dep(2)-ceil(L(2)/2)+L(2),1);
ZRtest = im(P_dep(1)-ceil(L(1)/2)+1-dP_ini_limit(1,1):P_dep(1)-ceil(L(1)/2)+L(1)+dP_ini_limit(1,2),...
    P_dep(2)-ceil(L(2)/2)+1-dP_ini_limit(2,1):P_dep(2)-ceil(L(2)/2)+L(2)+dP_ini_limit(2,2),2);

[pos_best_corr, vcorr]=CostFunctions(B1,ZRtest,[1 1],1);

dP_ini_limit=(pos_best_corr-1-dP_ini_limit(:,1)');


%progress bar
progress_bar=timerbar(0,'calculating ...');

%triangle scan of the image
for cercle_it=0:nb_boucle
    borne_min_N=max(1-N_mat_dep,-cercle_it);
    borne_max_N=min(N_max-N_mat_dep,cercle_it);
    %progress bar
    timerbar(nb_pts_calcul/nb_pts_est,progress_bar);
    nb_pts_calcul=nb_pts_calcul+(borne_max_N-borne_min_N+1);
    %triangle scan of the image
    for jN=borne_min_N:borne_max_N
        
        jM=min(cercle_it-abs(jN),M_max);
        clc
        
        for n = 1:(N-1) %images loop
            
            %translations initialisation
            if n==1 %estimation between the first two images
                
                if jM==0 & jN==0
                    
                    dP=dP_ini_limit;
                    
                else
                    
                    if jN~=0
                        dP(1)=round(champ_dep_est(jM+1,N_mat_dep+jN-sign(jN),4,1,n));
                    else
                        dP(1)=round(champ_dep_est(jM,N_mat_dep,4,1,n));
                    end
                    
                    if jM~=0
                        dP(2)=round(champ_dep_est(jM,N_mat_dep+jN,4,2,n));
                    else
                        dP(2)=round(champ_dep_est(1,N_mat_dep+jN-sign(jN),4,2,n));
                    end
                    
                end
            else %only for sequences...
                
                % the algorithme only estimates for image couples...
                
            end
            
            %Central point of the study zone
            if n == 1
                P=P_dep+Grid.*[jM jN];
            else
                % the algorithme only estimates for image couples...
            end
            
            %Current study zone on the first image
            if n==1
                Pat_ini=im(P(1)-L(1):P(1)+L(1),P(2)-L(2):P(2)+L(2),n);
            else
                % the algorithme only estimates for image couples...
            end
            
            
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            dP = [0 0];
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            
            ZR = im(P(1)-L(1)+dP(1):P(1)+L(1)+dP(1),P(2)-L(2)+dP(2):P(2)+L(2)+dP(2),n+1);
            
            phi1_ini = phi1(P(1)-L(1):P(1)+L(1),P(2)-L(2):P(2)+L(2));
            theta1_ini = theta1(P(1)-L(1):P(1)+L(1),P(2)-L(2):P(2)+L(2));
            freq1_ini = freq1(P(1)-L(1):P(1)+L(1),P(2)-L(2):P(2)+L(2));
            
            phi2_ZR = phi2(P(1)-L(1)+dP(1):P(1)+L(1)+dP(1),P(2)-L(2)+dP(2):P(2)+L(2)+dP(2));
            theta2_ZR = theta2(P(1)-L(1)+dP(1):P(1)+L(1)+dP(1),P(2)-L(2)+dP(2):P(2)+L(2)+dP(2));
            freq2_ZR = freq2(P(1)-L(1)+dP(1):P(1)+L(1)+dP(1),P(2)-L(2)+dP(2):P(2)+L(2)+dP(2));
            
            %estimation du champ de déplacement en P
            dep=bilinearEstimation...
                (Pat_ini,ZR,L,tR,nb_it,vI,t_cor,type_calcul_dep,...
                phi1_ini,theta1_ini,freq1_ini,phi2_ZR,theta2_ZR,freq2_ZR);
            
            %add the initial translations
            dep(4,:)=dep(4,:)+dP;
            %Bilinear parameters structure
            champ_dep_est(jM+1,N_mat_dep+jN,:,:,n)=dep;
            
        end
        
    end
end
close(progress_bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          ESTIMATION LOCALE         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chp_dep=bilinearEstimation...
    (PA,ZR,L,tR,nb_it,vI,type_correl,type_calcul_dep,phi1_ini,theta1_ini,freq1_ini,phi2_ZR,theta2_ZR,freq2_ZR)


% fonction qui estime localement les paramètres du modèle bilinéaire

%%% INPUT
% PA= PA = study zone on the initial image
% ZR= ZR = search zone on the second image
% L
% tR
% nb_it
% vI
% type_correl
% type_calcul_dep


% US_IRM

%%% OUTPUT
% chp_dep = estimated bilinear parameters for the current study zone


%%%%%%%%%%%% CONSTANTES %%%%%%%%%%%%%%%
%matrice calcul des déplacements
M_dep = 0.5* [-1/(L(2)-1)             1/(L(2)-1)             1/(L(2)-1)            -1/(L(2)-1) ; ...
    -1/(L(1)-1)            -1/(L(1)-1)             1/(L(1)-1)             1/(L(1)-1) ; ...
    2/((L(1)-1)*(L(2)-1)) -2/((L(1)-1)*(L(2)-1))  2/((L(1)-1)*(L(2)-1)) -2/((L(1)-1)*(L(2)-1)) ; ...
    1/2                    1/2                    1/2                    1/2               ];

%grids used to deform the study zone from one iteration to another
[X Y]=meshgrid(-L(2):L(2),-L(1):L(1));
XY=X.*Y;

%initial study zone
Pat_ini=PA;

%%%%%% Blocks translations estimation %%%%%%%%

for i=1:nb_it
    
    %matrice niveau de gris des patterns de l'image initiale + 1 pixel autour
    Pat_im1{1}=Pat_ini(2:L(1)+1,2:L(2)+1);
    Pat_im1{2}=Pat_ini(2:L(1)+1,L(2)+1:2*L(2));
    Pat_im1{3}=Pat_ini(L(1)+1:2*L(1),L(2)+1:2*L(2));
    Pat_im1{4}=Pat_ini(L(1)+1:2*L(1),2:L(2)+1);
    
    %calcul du déplacement sur chaque pattern
    try
        if type_correl == 6
            if i==1 % at the first iteration
                
                
                phi_im1{1}=phi1_ini(2:L(1)+1,2:L(2)+1);
                theta_im1{1}=theta1_ini(2:L(1)+1,2:L(2)+1);
                freq_im1{1}=freq1_ini(2:L(1)+1,2:L(2)+1);
                phi_im2{1}=phi2_ZR(2:L(1)+1,2:L(2)+1);
                theta_im2{1}=theta2_ZR(2:L(1)+1,2:L(2)+1);
                freq_im2{1}=freq2_ZR(2:L(1)+1,2:L(2)+1);
                
                phi_im1{2}=phi1_ini(2:L(1)+1,L(2)+1:2*L(2));
                theta_im1{2}=theta1_ini(2:L(1)+1,L(2)+1:2*L(2));
                freq_im1{2}=freq1_ini(2:L(1)+1,L(2)+1:2*L(2));
                phi_im2{2}=phi2_ZR(2:L(1)+1,L(2)+1:2*L(2));
                theta_im2{2}=theta2_ZR(2:L(1)+1,L(2)+1:2*L(2));
                freq_im2{2}=freq2_ZR(2:L(1)+1,L(2)+1:2*L(2));
                
                phi_im1{3}=phi1_ini(L(1)+1:2*L(1),L(2)+1:2*L(2));
                theta_im1{3}=theta1_ini(L(1)+1:2*L(1),L(2)+1:2*L(2));
                freq_im1{3}=freq1_ini(L(1)+1:2*L(1),L(2)+1:2*L(2));
                phi_im2{3}=phi2_ZR(L(1)+1:2*L(1),L(2)+1:2*L(2));
                theta_im2{3}=theta2_ZR(L(1)+1:2*L(1),L(2)+1:2*L(2));
                freq_im2{3}=freq2_ZR(L(1)+1:2*L(1),L(2)+1:2*L(2));
                
                phi_im1{4}=phi1_ini(L(1)+1:2*L(1),2:L(2)+1);
                theta_im1{4}=theta1_ini(L(1)+1:2*L(1),2:L(2)+1);
                freq_im1{4}=freq1_ini(L(1)+1:2*L(1),2:L(2)+1);
                phi_im2{4}=phi2_ZR(L(1)+1:2*L(1),2:L(2)+1);
                theta_im2{4}=theta2_ZR(L(1)+1:2*L(1),2:L(2)+1);
                freq_im2{4}=freq2_ZR(L(1)+1:2*L(1),2:L(2)+1);
                
                
                for patt=1:4
                    
                    [dPx dPy] = MonogenicAdjust...
                        (phi_im1{patt},theta_im1{patt},freq_im2{patt},phi_im2{patt},theta_im2{patt},freq_im2{patt});
                    dP{patt} = [dPy dPx];
                                        
                end
            else % after the first iteration
                Pat_im2{1}=ZR(2:L(1)+1,2:L(2)+1);
                Pat_im2{2}=ZR(2:L(1)+1,L(2)+1:2*L(2));
                Pat_im2{3}=ZR(L(1)+1:2*L(1),L(2)+1:2*L(2));
                Pat_im2{4}=ZR(L(1)+1:2*L(1),2:L(2)+1);
                
                for patt=1:4
                    
                    [dPx dPy] = MonogenicAdjust2(Pat_im1{patt},Pat_im2{patt});
                    dP{patt} = [dPy dPx];
                    
                end
            end
        else
            %             for patt=1:4
            %                 %grille d'interpolation
            %                 ixi_new=ixi+SP_ZR{patt}(2)-1;
            %                 iyi_new=iyi+SP_ZR{patt}(1)-1;
            %                 %Search zone interpolation
            %                 Pat_ZR=interp2(ix,iy,ZR,ixi_new,iyi_new);
            %
            %                 [pos_best_corr, vcorr{patt}]=CostFunctions(Pat_im1{patt},Pat_ZR,vI.^i,type_correl);
            %
            %                 dP{patt}=((pos_best_corr-1)./vI-tR(:,1)')./(vI.^(i-1))+dP{patt};
            %             end
        end
    catch
        break
    end
    
    %Bilinear parameters calculation
    dep=[dP{1} ; dP{2} ; dP{3} ; dP{4}];     %Blocks translations
    U=M_dep*dep(:,2);                        %[au bu cu du]
    V=M_dep*dep(:,1);                        %[av bv cv dv]
    
    % patterns deformation using estimated bilinear parameters at the previous
    % iteration
    if i~=nb_it
        XP=X+U(1)*X+U(2)*Y+U(3)*XY;           %deformation on X
        YP=Y+V(1)*X+V(2)*Y+V(3)*XY;           %deformation on Y
        Pat_ini=griddata(XP,YP,PA,X,Y);
        Pat_ini(find(isnan(Pat_ini)))=0;
    end
    
end
chp_dep=[V U];









