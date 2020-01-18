function [Vx Vy] = PhasePlanAdjustAffine(P1_r,P2_r,P1_s,P2_s,fx,fy,type_fenetre)

%   PhasePlanAdjustAffine Estimation locale des parameteres du modele de 
%   de deformation affine via 2 estimations par moindres carrées des 3x2
%   paramètres à partir de la somme et de la différence des différence de
%   phase obtenues après avoir calculé les signaux analytiques sur les 
%   images entières
%   L'estimations par moindres carrées est faites après avoir ajusté et
%   seuiller les valeurs correspondant sauts de phases(ref méthode AAPE) 

%   [dx dy] = PhaseAdjustAffine(r,s,fx,fy) computes the 2D estimated local
%   affine parameters 3x2 using phase adjustement and  two least square
%   estimators
%   Inputs : r,s 2D signals
%            fx,fy frequencies along each direction x and y
%   Outputs : Vx,Vy estimated local affine parameters "dx, ax, bx, dy, ay
%   and by"

%   Author: Florian Douziech, Mai 20011.


% Reference
%   [1] F.Douziech, A. Basarab, H. Liebgott Affine phase based motion 
%   estimation applied to echocardiography, IEEE Ultrasonics Symposium,
%   2011.

% Exemple

% x = 1:30; y = 1:30; [X,Y] = meshgrid(x,y);
% fx = 1/10; fy = 1/10;
% % True shifts
% dx = 0.45; dy = 0.45;
% 
% r = cos(2*pi*fx*X) .* cos(2*pi*fy*Y);
% s = cos(2*pi*fx*(X-dx)) .* cos(2*pi*fy*(Y-dy));
% [r1,r2] = ashahn(r); [s1,s2] = ashahn(s);
% P1_r=atan2(imag(r1),real(r1));
% P2_r=atan2(imag(r2),real(r2));
% P1_s=atan2(imag(s1),real(s1));
% P2_s=atan2(imag(s2),real(s2));
% figure; imagesc(x,y,r); title('r')
% figure; imagesc(x,y,s); title('s')
%
% [Vx_estim Vy_estim] = PhaseAdjustAffine(P1_r,P2_r,P1_s,P2_s,fx,fy)

%% TRAITEMENT DES DONNEES(AJUSTEMENT ET SEUILLAGE)
%% estimation grossiere des parametre dans une petite fenetre
%suffisament petite pour evité tout probleme de deroulement

tbx=5;      %taille de la ROI choisi pour ajuster les sauts de phase
tby=5;



        [y_Bmax x_Bmax] =size(P1_r);                  %parametre de taille du block dans lequel on estime le mouvement
                                                      %prendre des nombres impair pour des raison de symetrie
        center_y=round(y_Bmax/2);                   %centre de la ROI choisi au depart pour estimé le mouvement
        center_x=round(x_Bmax/2);                                               
                                                      
                                                      
        phi_diff1 = P1_r - P1_s;
        phi_diff2 = P2_r - P2_s;

        %%estimation parametrique(affine)
        seuil=pi;

        %systeme d'equation du modele affine: données 
        model_y=(phi_diff1+phi_diff2);
        model_x=(phi_diff1-phi_diff2);
        
       
        
        %systeme determinant le plan qui permet l'ajustement des saut de
        %phase
        model_y_Adjust=model_y(center_y-tby:center_y+tby,center_x-tbx:center_x+tbx);
        model_x_Adjust=model_x(center_y-tby:center_y+tby,center_x-tbx:center_x+tbx);
        %creation des grilles de coordonnées du bloc
        x_block_Ad=-tbx:tbx;
        y_block_Ad=-tby:tby;    
        [Xe Ye]=meshgrid(x_block_Ad,y_block_Ad);
        %seuillage pour eliminer les sauts de phase
        ind_saut1=(abs(model_y_Adjust)<seuil);       
        ind_saut2=(abs(model_x_Adjust)<seuil);  
        x_pos1=Xe(ind_saut1);
        y_pos1=Ye(ind_saut1); 
        x_pos2=Xe(ind_saut2);
        y_pos2=Ye(ind_saut2);
        model_ya=model_y_Adjust(ind_saut1);
        model_xa=model_x_Adjust(ind_saut2);

        D1=[ones(length(model_ya),1) x_pos1 y_pos1]; 
    
        Vad1=D1\model_ya;%vecteur des parametres   
        D2=[ones(length(model_xa),1) x_pos2 y_pos2]; 
        Vad2=D2\model_xa;%vecteur des parametres   
        
        %% reconstruction du plan qui sert de critere d'ajustement des sauts de phase
        f=ones((tbx*2+1)*(tbx*2+1),1);
        x_block=-(fix(x_Bmax/2)):(fix(x_Bmax/2));
        y_block=-(fix(y_Bmax/2)):(fix(y_Bmax/2)); 
        [Xe Ye]=meshgrid(x_block,y_block);

        C=[ones(x_Bmax*y_Bmax,1) reshape(Xe,x_Bmax*y_Bmax,1) reshape(Ye,x_Bmax*y_Bmax,1)];
        seuil_adjust1=C*Vad1;
        E1=reshape(seuil_adjust1,y_Bmax,x_Bmax);
        seuil_adjust2=C*Vad2;
        E2=reshape(seuil_adjust2,y_Bmax,x_Bmax);
        
        %% ajustement itératif suivant le nombre de plan à derouler
        %les donnees correspondant au 2 plans des parametres affines sont
        %decaler de +ou- 2*pi si elles exèdent l'intervalle [-2pi; 2pi]
        %l'opération et repeter tant qu'il reste des valeurs en dehors de
        %cet intervalle
        model_xs=model_x;
        while(sum(sum(model_xs>E2+2*pi))||sum(sum(model_xs<E2-2*pi)))
        seuilplus=(model_xs>E2+pi/2);
        seuilmoins=(model_xs<E2-pi/2);
        model_xs(seuilplus)=model_xs(seuilplus)-2*pi;
        model_xs(seuilmoins)=model_xs(seuilmoins)+2*pi;
        end
            clear seuilplus seuilmoins
        model_ys=model_y;
        while(sum(sum(model_ys>E1+2*pi))||sum(sum(model_ys<E1-2*pi)))
        seuilplus=(model_ys>E1+pi/2);
        seuilmoins=(model_ys<E1-pi/2);
        model_ys(seuilplus)=model_ys(seuilplus)-2*pi;
        model_ys(seuilmoins)=model_ys(seuilmoins)+2*pi;
        end
        clear seuilplus seuilmoins
        
        %% seuillage du plan pour eliminer les valeur aberrantes dues au
        %bruit
        seuilplus1=(model_ys>E1+pi/2)+(model_ys<E1-pi/2);
        seuilplus1=seuilplus1>0;


        model_ys(seuilplus1)=NaN;
        
        seuilplus2=(model_xs>E2+pi/2)+(model_xs<E2-pi/2);
        seuilplus2=seuilplus2>0;

        model_xs(seuilplus2)=NaN;

        
        
     %% ESTIMATION FINALE DES PARAMETRES APRES TRAITEMENT DES DONNEES   
        
        %fenetre d'apodisation 
        
        chaine1=[type_fenetre, '(x_Bmax)'];
        chaine2=[type_fenetre, '(y_Bmax)'];
%         fenetreHx=rectwin(x_Bmax);
%         fenetreHy=rectwin(y_Bmax);
        fenetreHx=eval(chaine1);
        fenetreHy=eval(chaine2);

        fenetreH=fenetreHy*fenetreHx';  
        %seuillage pour eliminer les sauts de phase
 
        %vecteurs de la matrice de regression pondéré
        %sans les points comportant des sauts de phase 
        x_pond=Xe.*fenetreH;
        y_pond=Ye.*fenetreH; 
        


        %pondération des données
        model1=model_ys.*fenetreH;
        model2=model_xs.*fenetreH;
        
     
        
        model_ys(isnan(model_ys))=[];
        [Lda Ldb]=size(model_ys);  
        Ld2=Lda*Ldb;
        donnee1=reshape(model_ys,Ld2,1);           
        seuilplus1=seuilplus1<1;        
        C3=[ones(Ld2,1) reshape(x_pond(seuilplus1),Ld2,1) reshape(y_pond(seuilplus1),Ld2,1)];
        Vy=C3\donnee1;
        
        model_xs(isnan(model_xs))=[];
        [Lda Ldb]=size(model_xs);  
        Ld2=Lda*Ldb;
        donnee2=reshape(model_xs,Ld2,1);           
        seuilplus2=seuilplus2<1;        
        C2=[ones(Ld2,1) reshape(x_pond(seuilplus2),Ld2,1) reshape(y_pond(seuilplus2),Ld2,1)];
        Vx=C2\donnee2;
        Vy = (1/(4*pi*fy)).*Vy;              
        Vx = (1/(4*pi*fx)).*Vx;   
       


  
        
end

                
     
                
          
           
  