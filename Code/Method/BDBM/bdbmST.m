function [champ_dep_est, x_grid, y_grid]=bdbmST(varargin);

% 2D parametric block matching motion estimation based on bilinear motion
% model

%INPUT
%  varargin{1}=im :    images im(:,:,:)
%  varargin{2}= Z_im :  estimated zone [ymin ymax xmin xmax]
%  varargin{3}=L :     ROIs size [Lv Lu]
%  varargin{4}=Grid :  mesh grid on the initial image [Gy Gx]
%  varargin{5}=vI :    interpolation factor for each iteration [int_y int_x]
%  varargin{6}=nb_it : number of iterations
%  varargin{7}=t_cor : type of local 2D spatial delay estimation.
%                                           1 = SAD ; 2 = normalized correlation; 3 = SAD + gradient descend search;
%                                           4 = elimination successive;
%                                           5 = phase least square estimation (LSPestimation) on the correlation function
%                                           6 = LSP on signals
%
%  varargin{8}=type_cacul_dep: facon de calculer le d�placement moyen lors du calcul des valeur estim�s. 1= moyennne,2 = m�diane
%  varargin{9}= dep_init : nodes displacements for images I1 In-1
%  varargin{10}=tR :   tolerance of search regions reported to the blocks size [yhaut ybas xgauche xdroite]
%  varargin{11}=f=[fy fx] : 2-D frequency vector
%  varargin{12}= US_IRM : 1 if RF, 0 if tagged RMI

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
%       Elsevier, Medical Image Analysis, 2008.



a=varargin{4};
if ~exist('a','var')
    varargin{4} = varargin{3};
end
a=varargin{5};
if ~exist('a','var')
    varargin{5} = [3 3];
end
a=varargin{6};
if ~exist('a','var')
    varargin{6} = 2;
end
a=varargin{7};
if ~exist('a','var')
    varargin{7} = 1;
end
a=varargin{10};
if ~exist('a','var')
    varargin{10}=[1  1 ; 1 1];
end

if (varargin{7} == 6) && (varargin{12} == 1)
    varargin{5} = [1 1]; % no interpolation in the case of LSE phase estimator
    varargin{10} = [1 1 ; 1 1];

    [R1,R2]=ashahn(varargin{1}(:,:,1));
    [S1,S2]=ashahn(varargin{1}(:,:,2));
    P1_im1=atan2(imag(R1),real(R1));
    P2_im1=atan2(imag(R2),real(R2));
    P1_im2=atan2(imag(S1),real(S1));
    P2_im2=atan2(imag(S2),real(S2));
end


%Numbers of images in the sequence
N = size(varargin{1},3);

global champ_dep_est
global x_grid
global y_grid

%Search zone size
tZ=varargin{3}-1+sum(varargin{10}');
%Study zone size
tROI=2*varargin{3}-2+sum(varargin{10}');

%taille de la zone d'�tude
long_Z_im=diff(varargin{2}(1,:))-tROI(1);
M_max=fix(long_Z_im/varargin{4}(1))+1;
larg_Z_im=diff(varargin{2}(2,:))-tROI(2);
N_max=fix(larg_Z_im/varargin{4}(2));
%point de calcul le plus � gauche
debut=tZ(2)-varargin{10}(2,2)+varargin{2}(2,1)+fix((larg_Z_im-(N_max-1)*varargin{4}(2))/2);
%centre de d�part
centre_im_N=fix(N_max/2)*varargin{4}(2)+debut;     %image coordinates
N_mat_dep=fix(N_max/2)+1;                   %coordinates in the matrix

%first node to be estimated
P_dep=[varargin{3}(1)+varargin{10}(1,1)+varargin{2}(1,1) centre_im_N];

%grid of the nodes to be estimated
x_grid=debut:varargin{4}(2):debut+(N_max-1)*varargin{4}(2);
y_grid=P_dep(1):varargin{4}(1):P_dep(1)+M_max*varargin{4}(1);

%number of nodes to be estimated
nb_pts_est=(N_max-1)*M_max;


nb_pts_calcul=0;

%number of loops
nb_boucle=N_max+M_max-N_mat_dep+1;

champ_dep_est=NaN*ones(M_max,N_max,4,2);

%progress bar
progress_bar=timerbar(0,'calculating ...');

%triangle scan of the image
for cercle_it=0:nb_boucle
    borne_min_N=max(1-N_mat_dep,-cercle_it);
    borne_max_N=min(N_max-N_mat_dep,cercle_it);
    %progress bar
%     timerbar(nb_pts_calcul/nb_pts_est,progress_bar);
    nb_pts_calcul=nb_pts_calcul+(borne_max_N-borne_min_N+1);
    %triangle scan of the image
    for jN=borne_min_N:borne_max_N

        jM=min(cercle_it-abs(jN),M_max);
        clc
        %translations initialisation
        dP(1) = round(varargin{9}(jM+1,N_mat_dep+jN,4,1));
        dP(2) = round(varargin{9}(jM+1,N_mat_dep+jN,4,2));


        %Central point of the study zone
        P=P_dep+varargin{4}.*[jM jN];
        Pat_ini=varargin{1}(P(1)-varargin{3}(1):P(1)+varargin{3}(1),P(2)-varargin{3}(2):P(2)+varargin{3}(2),1);


        if varargin{7} == 6 % LSP estimator directly on the signals
            ZR = varargin{1}(P(1)-varargin{3}(1)+dP(1):P(1)+varargin{3}(1)+dP(1),P(2)-varargin{3}(2)+dP(2):P(2)+varargin{3}(2)+dP(2),2);
            if varargin{12}==1
                P1_ini = P1_im1(P(1)-varargin{3}(1):P(1)+varargin{3}(1),P(2)-varargin{3}(2):P(2)+varargin{3}(2));
                P2_ini = P2_im1(P(1)-varargin{3}(1):P(1)+varargin{3}(1),P(2)-varargin{3}(2):P(2)+varargin{3}(2));
                P1_ZR = P1_im2(P(1)-varargin{3}(1)+dP(1):P(1)+varargin{3}(1)+dP(1),P(2)-varargin{3}(2)+dP(2):P(2)+varargin{3}(2)+dP(2));
                P2_ZR = P2_im2(P(1)-varargin{3}(1)+dP(1):P(1)+varargin{3}(1)+dP(1),P(2)-varargin{3}(2)+dP(2):P(2)+varargin{3}(2)+dP(2));
            else
                % for IRM we calculate the phase locally, as we have to
                % rescale the signals around 0
                P1_ini = 0; P2_ini = 0; P1_ZR = 0; P2_ZR = 0;
            end
            %estimation du champ de d�placement en P
            dep=bilinearEstimation(Pat_ini,ZR,varargin{3},varargin{10},varargin{6},varargin{5},varargin{7},varargin{8},varargin{11},P1_ini,P2_ini,P1_ZR,P2_ZR,varargin{12});
        else
            %Search zone
            PsZ=P-varargin{3}+1-varargin{10}(:,1)'+dP;     %Up left pixel of the search zone
            ZR=varargin{1}(PsZ(1):PsZ(1)+tROI(1),PsZ(2):PsZ(2)+tROI(2),2);  %Search zone on the second image
            %estimation du champ de d�placement en P
            dep=bilinearEstimation(Pat_ini,ZR,varargin{3},varargin{10},varargin{6},varargin{5},varargin{7},varargin{8},varargin{11});
        end

        %             %estimation du champ de d�placement en P
        %             dep=bilinearEstimation(Pat_ini,ZR,varargin{3},varargin{10},varargin{6},varargin{5},varargin{7},varargin{8},varargin{11});
        %add the initial translations
        dep(4,:)=dep(4,:)+dP;
        %Bilinear parameters structure
        champ_dep_est(jM+1,N_mat_dep+jN,:,:)=dep;


    end
end
close(progress_bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          ESTIMATION LOCALE         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chp_dep=bilinearEstimation(varargin)

% fonction qui estime localement les param�tres du mod�le bilin�aire

%%% INPUT
% varargin{1}= PA = study zone on the initial image
% varargin{2}= ZR = search zone on the second image
% varargin{3}= L
% varargin{4}= tR
% varargin{5}= nb_it
% varargin{6}= vI
% varargin{7}= type_correl
% varargin{8}= type_calcul_dep
% varargin{9}= [fy fx] = 2-D frequency vector
% PLANS pour estimation LSPs at the first iteration
% varargin{10}= P1_ini
% varargin{11}= P2_ini
% varargin{12}= P1_ZR
% varargin{13}= P2_ZR
% varargin{14} = US_IRM


%%% OUTPUT
% chp_dep = estimated bilinear parameters for the current study zone


%%%%%%%%%%%% CONSTANTES %%%%%%%%%%%%%%%
%matrice calcul des d�placements
M_dep = 0.5* [-1/(varargin{3}(2)-1)             1/(varargin{3}(2)-1)             1/(varargin{3}(2)-1)            -1/(varargin{3}(2)-1) ; ...
    -1/(varargin{3}(1)-1)            -1/(varargin{3}(1)-1)             1/(varargin{3}(1)-1)             1/(varargin{3}(1)-1) ; ...
    2/((varargin{3}(1)-1)*(varargin{3}(2)-1)) -2/((varargin{3}(1)-1)*(varargin{3}(2)-1))  2/((varargin{3}(1)-1)*(varargin{3}(2)-1)) -2/((varargin{3}(1)-1)*(varargin{3}(2)-1)) ; ...
    1/2                    1/2                    1/2                    1/2               ];

%grids used to deform the study zone from one iteration to another
[X Y]=meshgrid(-varargin{3}(2):varargin{3}(2),-varargin{3}(1):varargin{3}(1));
XY=X.*Y;

%initial study zone
Pat_ini=varargin{1};

if varargin{7} == 5 | varargin{7} == 6
    % no interpolation with LSE phase estimator
else
    %Search zone size
    tZ=varargin{3}-1+sum(varargin{4}');
    sumtR=sum(varargin{4}');

    % initial translations
    dP{1}=[0 0]; dP{2}=[0 0]; dP{3}=[0 0]; dP{4}=[0 0];
    dep=[dP{1}; dP{2}; dP{3} ;dP{4}];

    %up left corners of each of the 4 blocks of the study zone
    SP_ZR{1}=[1 1];
    SP_ZR{2}=[1 varargin{3}(2)];
    SP_ZR{3}=[varargin{3}(1) varargin{3}(2)];
    SP_ZR{4}=[varargin{3}(1) 1];
    %interpolation grid of the search region
    [ix,iy]=meshgrid(1:varargin{3}(2)+tZ(2),1:varargin{3}(1)+tZ(1));
end

%%%%%% Blocks translations estimation %%%%%%%%

for i=1:varargin{5}

    %matrice niveau de gris des patterns de l'image initiale + 1 pixel autour
    Pat_im1{1}=Pat_ini(2:varargin{3}(1)+1,2:varargin{3}(2)+1);
    Pat_im1{2}=Pat_ini(2:varargin{3}(1)+1,varargin{3}(2)+1:2*varargin{3}(2));
    Pat_im1{3}=Pat_ini(varargin{3}(1)+1:2*varargin{3}(1),varargin{3}(2)+1:2*varargin{3}(2));
    Pat_im1{4}=Pat_ini(varargin{3}(1)+1:2*varargin{3}(1),2:varargin{3}(2)+1);

    if varargin{7} == 6

    else
        % mean value of the estimated translations of the 4 blocks
        switch(varargin{8})
            case 1
                %en prenant la valeur du d�placment moyen
                depm=mean(dep);
            case 2
                %en prenant la valeur m�diane du d�placement
                depm=median(dep);
        end
        pZR=1+depm+(1-1./(varargin{6}.^(i-1))).*varargin{4}(:,1)';
        %grille d'interpolation de la nouvelle zone de recherche
        [ixi,iyi]=meshgrid(pZR(2):1/(varargin{6}(2)^i):pZR(2)+varargin{3}(2)-1+sumtR(2)/(varargin{6}(2)^(i-1)),pZR(1):1/(varargin{6}(1)^i):pZR(1)+varargin{3}(1)-1+sumtR(1)/(varargin{6}(1)^(i-1)));
    end

    %calcul du d�placement sur chaque pattern
    try
        if varargin{7} == 6
            if i==1 % at the first iteration
                if varargin{14} == 1
                    P1_im1{1}=varargin{10}(2:varargin{3}(1)+1,2:varargin{3}(2)+1);
                    P2_im1{1}=varargin{11}(2:varargin{3}(1)+1,2:varargin{3}(2)+1);
                    P1_im2{1}=varargin{12}(2:varargin{3}(1)+1,2:varargin{3}(2)+1);
                    P2_im2{1}=varargin{13}(2:varargin{3}(1)+1,2:varargin{3}(2)+1);

                    P1_im1{2}=varargin{10}(2:varargin{3}(1)+1,varargin{3}(2)+1:2*varargin{3}(2));
                    P2_im1{2}=varargin{11}(2:varargin{3}(1)+1,varargin{3}(2)+1:2*varargin{3}(2));
                    P1_im2{2}=varargin{12}(2:varargin{3}(1)+1,varargin{3}(2)+1:2*varargin{3}(2));
                    P2_im2{2}=varargin{13}(2:varargin{3}(1)+1,varargin{3}(2)+1:2*varargin{3}(2));

                    P1_im1{3}=varargin{10}(varargin{3}(1)+1:2*varargin{3}(1),varargin{3}(2)+1:2*varargin{3}(2));
                    P2_im1{3}=varargin{11}(varargin{3}(1)+1:2*varargin{3}(1),varargin{3}(2)+1:2*varargin{3}(2));
                    P1_im2{3}=varargin{12}(varargin{3}(1)+1:2*varargin{3}(1),varargin{3}(2)+1:2*varargin{3}(2));
                    P2_im2{3}=varargin{13}(varargin{3}(1)+1:2*varargin{3}(1),varargin{3}(2)+1:2*varargin{3}(2));

                    P1_im1{4}=varargin{10}(varargin{3}(1)+1:2*varargin{3}(1),2:varargin{3}(2)+1);
                    P2_im1{4}=varargin{11}(varargin{3}(1)+1:2*varargin{3}(1),2:varargin{3}(2)+1);
                    P1_im2{4}=varargin{12}(varargin{3}(1)+1:2*varargin{3}(1),2:varargin{3}(2)+1);
                    P2_im2{4}=varargin{13}(varargin{3}(1)+1:2*varargin{3}(1),2:varargin{3}(2)+1);

                    for patt=1:4
                        [dPx dPy] = PhaseAdjust3(P1_im1{patt},P2_im1{patt},P1_im2{patt},P2_im2{patt},varargin{9}(2),varargin{9}(1));
                        dP{patt} = [dPy dPx];
                    end
                else
                    Pat_im2{1}=varargin{2}(2:varargin{3}(1)+1,2:varargin{3}(2)+1);
                    Pat_im2{2}=varargin{2}(2:varargin{3}(1)+1,varargin{3}(2)+1:2*varargin{3}(2));
                    Pat_im2{3}=varargin{2}(varargin{3}(1)+1:2*varargin{3}(1),varargin{3}(2)+1:2*varargin{3}(2));
                    Pat_im2{4}=varargin{2}(varargin{3}(1)+1:2*varargin{3}(1),2:varargin{3}(2)+1);

                    for patt=1:4
                       [dPx dPy] = PhaseAdjust2(Pat_im1{patt},Pat_im2{patt},varargin{9}(2),varargin{9}(1),varargin{14});
                        dP{patt} = [dPy dPx];

                    end
                end
            else % after the first iteration
                Pat_im2{1}=varargin{2}(2:varargin{3}(1)+1,2:varargin{3}(2)+1);
                Pat_im2{2}=varargin{2}(2:varargin{3}(1)+1,varargin{3}(2)+1:2*varargin{3}(2));
                Pat_im2{3}=varargin{2}(varargin{3}(1)+1:2*varargin{3}(1),varargin{3}(2)+1:2*varargin{3}(2));
                Pat_im2{4}=varargin{2}(varargin{3}(1)+1:2*varargin{3}(1),2:varargin{3}(2)+1);

                for patt=1:4

                    [dPx dPy] = PhaseAdjust2(Pat_im1{patt},Pat_im2{patt},varargin{9}(2),varargin{9}(1),varargin{14});
                    dP{patt} = [dPy dPx];

                end
            end
        else
            for patt=1:4
                %grille d'interpolation
                ixi_new=ixi+SP_ZR{patt}(2)-1;
                iyi_new=iyi+SP_ZR{patt}(1)-1;
                %Search zone interpolation
                Pat_ZR=interp2(ix,iy,varargin{2},ixi_new,iyi_new);

                [pos_best_corr, vcorr{patt}]=CostFunctions(Pat_im1{patt},Pat_ZR,varargin{6}.^i,varargin{7});

                dP{patt}=((pos_best_corr-1)./varargin{6}-varargin{4}(:,1)')./(varargin{6}.^(i-1))+dP{patt};
            end
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
    if i~=varargin{5}
        XP=X+U(1)*X+U(2)*Y+U(3)*XY;           %deformation on X
        YP=Y+V(1)*X+V(2)*Y+V(3)*XY;           %deformation on Y
        Pat_ini=griddata(XP,YP,varargin{1},X,Y);
        Pat_ini(find(isnan(Pat_ini)))=0;
    end

end
chp_dep=[V U];

