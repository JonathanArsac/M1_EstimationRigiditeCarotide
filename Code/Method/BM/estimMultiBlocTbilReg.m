clear all
close all
clc

cd('/home/foxiler/git/M1_EstimationRigiditeCarotide/Donnees/Res_img/')

x = [66:75;112:121;68:77]; y = [100:109;74:83;65:74]; %path2

ROI(1) = size(y,2); ROI(2) = size(x,2);

M = 0.5* [-1/(ROI(2)-1) 1/(ROI(2)-1) 1/(ROI(2)-1) -1/(ROI(2)-1) ; ...
    -1/(ROI(1)-1) -1/(ROI(1)-1) 1/(ROI(1)-1) 1/(ROI(1)-1) ; ...
    2/((ROI(1)-1)*(ROI(2)-1)) -2/((ROI(1)-1)*(ROI(2)-1))  2/((ROI(1)-1)*(ROI(2)-1)) -2/((ROI(1)-1)*(ROI(2)-1)) ; ...
    1/2                    1/2                    1/2                    1/2               ];

%%%% Estimation using local rigid modtion or BDBM
% 1 if BM, 0 if BDBM
BM_BDBM = 0;

%%%%%%parametres PBM%%%%%%%%
Z_im=[ 20  150 ;   20   150];
Grid=[3 3]; L=[6 6]; vI=[1 1];
dP_ini_limit = [0 0 ; 0 0]; tR=[1 1 ;1 1];
type_estim_local = 6;
%%%%%%parametres BDBM%%%%%%%%
type_calcul_dep= 1; type_interp=1; nb_it = 1;


for k = 0:34

    %%%% charger les images et le vrai deplacement
    if k<10
        filedispx = ['disp_x_000'];
        filedispy = ['disp_y_000'];
        fileim = ['seq000'];
    else
        filedispx = ['disp_x_00'];
        filedispy = ['disp_y_00'];
        fileim = ['seq00'];
    end
    filename = [filedispx,num2str(k),'.img'];
    dispx = readimg(filename,160,160);
    filename = [filedispy,num2str(k),'.img'];
    dispy = readimg(filename,160,160);
    filename = [fileim,num2str(k),'.img'];
    a = readimg(filename,160,160);
    u(:,:,k+1) = dispx; v(:,:,k+1) = dispy; im(:,:,k+1) = a;

    %% d�formation d'un bloc au cours du temps (vraie)

    imagesc(im(:,:,k+1)); colormap gray;
    for B=1:size(x,1)
        d1_v(k+1,1,B) = dispx(y(B,1),x(B,1));d1_v(k+1,2,B) = dispx(y(B,1),x(B,end));d1_v(k+1,3,B) = dispx(y(B,end),x(B,end));d1_v(k+1,4,B) = dispx(y(B,end),x(B,1));
        d2_v(k+1,1,B) = dispy(y(B,1),x(B,1));d2_v(k+1,2,B) = dispy(y(B,1),x(B,end));d2_v(k+1,3,B) = dispy(y(B,end),x(B,end));d2_v(k+1,4,B) = dispy(y(B,end),x(B,1));

        %% line(x1, x2, y1, y2)
        line([x(B,1)+d1_v(k+1,1,B) x(B,end)+d1_v(k+1,2,B)], [y(B,1)+d2_v(k+1,1,B) y(B,1)+d2_v(k+1,2,B)],'Color','w')
        line([x(B,end)+d1_v(k+1,2,B) x(B,end)+d1_v(k+1,3,B)], [y(B,1)+d2_v(k+1,2,B) y(B,end)+d2_v(k+1,3,B)],'Color','w')
        line([x(B,end)+d1_v(k+1,3,B) x(B,1)+d1_v(k+1,4,B)], [y(B,end)+d2_v(k+1,3,B) y(B,end)+d2_v(k+1,4,B)],'Color','w')
        line([x(B,1)+d1_v(k+1,4,B) x(B,1)+d1_v(k+1,1,B)], [y(B,end)+d2_v(k+1,4,B) y(B,1)+d2_v(k+1,1,B)],'Color','w')      

        for p=1:4
            dP_v{p} = [d2_v(k+1,p,B) d1_v(k+1,p,B)];
        end
        dep_v=[dP_v{1} ; dP_v{2} ; dP_v{3} ; dP_v{4}];     %Blocks translations
        U_v(:,k+1,B)=M*dep_v(:,2);                        %[au bu cu du]
        V_v(:,k+1,B)=M*dep_v(:,1);
    end
    G(k+1)=getframe


    %% d�formation d'un bloc au cours du temps (estim�e)

    if k == 0
        % avec initialisation spatiale
        if BM_BDBM
            [champ_dep_est, x_grid, y_grid]=BM(im(:,:,1),im(:,:,k+1),Z_im,L,Grid,vI,type_estim_local,[0 0;0 0],tR,[0.1 0.1]);
        else
            Images(:,:,1) = im(:,:,1); Images(:,:,2) = im(:,:,k+1);
            [champ_dep_est x_grid y_grid] = bdbm(Images,Z_im,L,Grid,vI,nb_it,type_estim_local,type_calcul_dep,[0 0;0 0],tR,[0.1 0.1]);
        end
    else
        % avec initialisation temporelle
        if BM_BDBM
            [champ_dep_est, x_grid, y_grid]=BMST(im(:,:,1),im(:,:,k+1),champ_dep_est,Z_im,L,Grid,vI,type_estim_local,tR,[0.1 0.1],0);
        else
            Images(:,:,1) = im(:,:,1); Images(:,:,2) = im(:,:,k+1);
            [champ_dep_est, x_grid, y_grid] = bdbmST(Images,Z_im,L,Grid,vI,nb_it,type_estim_local,type_calcul_dep,champ_dep_est,tR,[0.1 0.1],0);
        end
    end    
    
    %%% champ dense
    if BM_BDBM
        depl_ax = champ_dep_est(:,:,1);
        depl_lat = champ_dep_est(:,:,2);
        [X_grid Y_grid] = meshgrid(x_grid,y_grid);
        [pts_lat pts_ax] = meshgrid(x_grid(1):x_grid(end),y_grid(1):y_grid(end));
        chp_dep_int(:,:,1) = interp2(X_grid,Y_grid,depl_ax,pts_lat,pts_ax,'cubic');
        chp_dep_int(:,:,2) = interp2(X_grid,Y_grid,depl_lat,pts_lat,pts_ax,'cubic');
    else
        [chp_dep_int,pts_ax,pts_lat] = denseField(champ_dep_est,x_grid, y_grid,Grid,type_interp);
    end
    
    imagesc(im(:,:,k+1)); colormap gray;

    for B=1:size(x,1) % pour chaque bloc affich�

        d1(k+1,1,B)=mean2(chp_dep_int((y(B,1)-L(2)/2-y_grid(1)):(y(B,1)+L(2)/2-y_grid(1)),(x(B,1)-L(1)/2-x_grid(1)):(x(B,1)+L(1)/2-x_grid(1)),2));
        d1(k+1,2,B)=mean2(chp_dep_int((y(B,1)-L(2)/2-y_grid(1)):(y(B,1)+L(2)/2-y_grid(1)),(x(B,end)-L(1)/2-x_grid(1)):(x(B,end)+L(1)/2-x_grid(1)),2));
        d1(k+1,3,B)=mean2(chp_dep_int((y(B,end)-L(2)/2-y_grid(1)):(y(B,end)+L(2)/2-y_grid(1)),(x(B,end)-L(1)/2-x_grid(1)):(x(B,end)+L(1)/2-x_grid(1)),2));
        d1(k+1,4,B)=mean2(chp_dep_int((y(B,end)-L(2)/2-y_grid(1)):(y(B,end)+L(2)/2-y_grid(1)),(x(B,1)-L(1)/2-x_grid(1)):(x(B,1)+L(1)/2-x_grid(1)),2));

        d2(k+1,1,B)=mean2(chp_dep_int((y(B,1)-L(2)/2-y_grid(1)):(y(B,1)+L(2)/2-y_grid(1)),(x(B,1)-L(1)/2-x_grid(1)):(x(B,1)+L(1)/2-x_grid(1)),1));
        d2(k+1,2,B)=mean2(chp_dep_int((y(B,1)-L(2)/2-y_grid(1)):(y(B,1)+L(2)/2-y_grid(1)),(x(B,end)-L(1)/2-x_grid(1)):(x(B,end)+L(1)/2-x_grid(1)),1));
        d2(k+1,3,B)=mean2(chp_dep_int((y(B,end)-L(2)/2-y_grid(1)):(y(B,end)+L(2)/2-y_grid(1)),(x(B,end)-L(1)/2-x_grid(1)):(x(B,end)+L(1)/2-x_grid(1)),1));
        d2(k+1,4,B)=mean2(chp_dep_int((y(B,end)-L(2)/2-y_grid(1)):(y(B,end)+L(2)/2-y_grid(1)),(x(B,1)-L(1)/2-x_grid(1)):(x(B,1)+L(1)/2-x_grid(1)),1));

        for p=1:4
            dP{p} = [d2(k+1,p,B) d1(k+1,p,B)];
        end

        dep=[dP{1} ; dP{2} ; dP{3} ; dP{4}];     %Blocks translations
        U(:,k+1,B)=M*dep(:,2);                        %[au bu cu du]
        V(:,k+1,B)=M*dep(:,1);
        line([x(B,1)+d1(k+1,1,B) x(B,end)+d1(k+1,2,B)], [y(B,1)+d2(k+1,1,B) y(B,1)+d2(k+1,2,B)],'Color','w')
        line([x(B,end)+d1(k+1,2,B) x(B,end)+d1(k+1,3,B)], [y(B,1)+d2(k+1,2,B) y(B,end)+d2(k+1,3,B)],'Color','w')
        line([x(B,end)+d1(k+1,3,B) x(B,1)+d1(k+1,4,B)], [y(B,end)+d2(k+1,3,B) y(B,end)+d2(k+1,4,B)],'Color','w')
        line([x(B,1)+d1(k+1,4,B) x(B,1)+d1(k+1,1,B)], [y(B,end)+d2(k+1,4,B) y(B,1)+d2(k+1,1,B)],'Color','w')
        
        rectangle('Position',[x(B,1)-L(1)/2+d1(k+1,1,B),y(B,1)-L(2)/2+d2(k+1,1,B),L(1),L(2)],'EdgeColor','r')
        rectangle('Position',[x(B,end)-L(1)/2+d1(k+1,2,B),y(B,1)-L(2)/2+d2(k+1,2,B),L(1),L(2)],'EdgeColor','r')
        rectangle('Position',[x(B,end)-L(1)/2+d1(k+1,3,B),y(B,end)-L(2)/2+d2(k+1,3,B),L(1),L(2)],'EdgeColor','r')
        rectangle('Position',[x(B,1)-L(1)/2+d1(k+1,4,B),y(B,end)-L(2)/2+d2(k+1,4,B),L(1),L(2)],'EdgeColor','r')

    end
    H(k+1)=getframe
end

t=0:(size(d1,1)-1);

for B=1:size(x,1)
    figure; plot3(x(B,1)+d1(:,1,B),y(B,1)+d2(:,1,B),t,'*-');
    hold;
    plot3(x(B,end)+d1(:,2,B),y(B,1)+d2(:,2,B),t,'r*-');
    plot3(x(B,end)+d1(:,3,B),y(B,end)+d2(:,3,B),t,'k*-');
    plot3(x(B,1)+d1(:,4,B),y(B,end)+d2(:,4,B),t,'m*-');

    plot3(x(B,1)+d1_v(:,1,B),y(B,1)+d2_v(:,1,B),t,'x--');
    plot3(x(B,end)+d1_v(:,2,B),y(B,1)+d2_v(:,2,B),t,'rx--');
    plot3(x(B,end)+d1_v(:,3,B),y(B,end)+d2_v(:,3,B),t,'kx--');
    plot3(x(B,1)+d1_v(:,4,B),y(B,end)+d2_v(:,4,B),t,'mx--');
end

% figure; movie(G,10,10)
% figure; movie(H,10,10)

