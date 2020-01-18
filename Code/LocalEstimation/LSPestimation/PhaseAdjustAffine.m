function [Vx Vy] = PhaseAdjustAffine(P1_r,P2_r,P1_s,P2_s,fx,fy,type_fenetre)

%   PhaseAdjustAffine  Estimation locale des parameteres du modele de 
%   de deformation affine via 2 estimations par moindres carrées des 3x2
%   paramètres à partir de la somme et de la différence des différence de
%   phase obtenues après avoir calculé les signaux analytiques sur les 
%   images entières
%   les données de 'model' en dehors de l'intervalle [-pi;+pi]
%   sont considéré comme fausse et sont éliminées de l'estimation

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
        [y_Bmax x_Bmax] =size(P1_r);                  %parametre de taille du block dans lequel on estime le mouvement
                                                      %prendre des nombres impair pour des raison de symetrie

        phi_diff1 = P1_r - P1_s;
        phi_diff2 = P2_r - P2_s;

        %%estimation parametrique(affine)
        seuil=pi;

        %systeme d'equation du modele affine: données 
        model_y=(phi_diff1+phi_diff2);
        model_x=(phi_diff1-phi_diff2);
        %fenetre d'apodisation 
        
        chaine1=[type_fenetre, '(x_Bmax)'];
        chaine2=[type_fenetre, '(y_Bmax)'];
%         fenetreHx=rectwin(x_Bmax);
%         fenetreHy=rectwin(y_Bmax);
        fenetreHx=eval(chaine1);
        fenetreHy=eval(chaine2);

        fenetreH=fenetreHy*fenetreHx';
        %creation des grilles de coordonnées du bloc
        x_block=-(fix(x_Bmax/2)):(fix(x_Bmax/2));
        y_block=-(fix(y_Bmax/2)):(fix(y_Bmax/2));        
        [Xe Ye]=meshgrid(x_block,y_block);
        %seuillage pour eliminer les sauts de phase
        ind_ssaut1=(abs(model_y)<seuil);       
        ind_ssaut2=(abs(model_x)<seuil);  
        %vecteurs de la matrice de regression pondéré
        %sans les points comportant des sauts de phase
        fenetre1=fenetreH(ind_ssaut1);  
        x_block1=Xe(ind_ssaut1).*fenetre1;
        y_block1=Ye(ind_ssaut1).*fenetre1; 
        
        fenetre2=fenetreH(ind_ssaut2);   
        x_block2=Xe(ind_ssaut2).*fenetre2;
        y_block2=Ye(ind_ssaut2).*fenetre2;

        %pondération des données
        model1=model_y(ind_ssaut1).*fenetre1;
        model2=model_x(ind_ssaut2).*fenetre2;
        
%% Least square estimator
        B1=[fenetre1 x_block1 y_block1];  
        BB1=inv(B1'*B1);
        Vy = BB1 *B1'*model1;%vecteur des parametres

        B2=[fenetre2 x_block2 y_block2];
        BB2=inv(B2'*B2);
        Vx = BB2 *B2'*model2;%vecteur des parametres modulo (1/(4*pi*fy))
            
        
%parametre estimé dans le bloc couran
        Vy = (1/(4*pi*fy)).*Vy;              
        Vx = (1/(4*pi*fx)).*Vx;        


  
        
end

                
     
                
          
           
  