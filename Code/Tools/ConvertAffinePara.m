function [chp_dep_int]=ConvertAffinePara(champ_dep_est,x_grid,y_grid,L,Tim,pixel_size_ax,pixel_size_lat)


%   ConvertAffinePara Converti les parametres du modele affine estimés dans 
%   chaque bloc en une carte de mouvement à la taille de l'image considérée


%   [chp_dep_int]=ConvertAffinePara(champ_dep_est,x_grid,y_grid,L,Tim,pixel_size_ax,pixel_size_lat)
%   convert estimated affine parameters into motion map 
%
%   Inputs : champ_dep_est  estimated affine parameters for each block 
%   x_grid,y_grid estimated block center in image
%   L  block size
%   Tim  image size
%   pixel_size_ax,pixel_size_lat  pixel size in mm
%
%   Outputs : chp_dep_int   motion map

%   Author: Florian Douziech, Mai 2011.




%taille des carte de parametre
[H_max K_max]=size(champ_dep_est(:,:,1,1));
%initialisation de la carte de champ dense de mouvement
%(de la tialle del'image)
chp_dep_int=zeros(Tim(1),Tim(2),3);


for jH=1:H_max
    for jK=1:K_max
            
        %grille des coordonnées d'un bloc
            x_block=-(fix(L(2)/2)):(fix(L(2)/2));
            y_block=-(fix(L(1)/2)):(fix(L(1)/2));
            [Xb_grid Yb_grid]=meshgrid(x_block,y_block);

    
        %matrice de regression(vecteur    
            B=[ones(L(1)*L(2),1) reshape(Xb_grid,L(1)*L(2),1) reshape(Yb_grid,L(1)*L(2),1)];

            
            dep_x_0 = B*squeeze(champ_dep_est(jH,jK,:,2));  
            dep_y_0 = B*squeeze(champ_dep_est(jH,jK,:,1));          

            dep_x = reshape(dep_x_0,L(2),L(1));        
            dep_y = reshape(dep_y_0,L(2),L(1));

            %sommation des mouvements de chaque pixel reconstruit à partir des parametres pour chaque blocs dans
            %dans la carte de champ dense de mouvement
            chp_dep_int(y_grid(jH)-ceil(L(1)/2)+1:y_grid(jH)-ceil(L(1)/2)+L(1),x_grid(jK)-ceil(L(2)/2)+1:x_grid(jK)-ceil(L(2)/2)+L(2),1)=...
                chp_dep_int(y_grid(jH)-ceil(L(1)/2)+1:y_grid(jH)-ceil(L(1)/2)+L(1),x_grid(jK)-ceil(L(2)/2)+1:x_grid(jK)-ceil(L(2)/2)+L(2),1)+dep_y';
            chp_dep_int(y_grid(jH)-ceil(L(1)/2)+1:y_grid(jH)-ceil(L(1)/2)+L(1),x_grid(jK)-ceil(L(2)/2)+1:x_grid(jK)-ceil(L(2)/2)+L(2),2)=...
                chp_dep_int(y_grid(jH)-ceil(L(1)/2)+1:y_grid(jH)-ceil(L(1)/2)+L(1),x_grid(jK)-ceil(L(2)/2)+1:x_grid(jK)-ceil(L(2)/2)+L(2),2)+dep_x';
            %nombre de dep ajouter aux pixels considérés
            chp_dep_int(y_grid(jH)-ceil(L(1)/2)+1:y_grid(jH)-ceil(L(1)/2)+L(1),x_grid(jK)-ceil(L(2)/2)+1:x_grid(jK)-ceil(L(2)/2)+L(2),3)=...
                chp_dep_int(y_grid(jH)-ceil(L(1)/2)+1:y_grid(jH)-ceil(L(1)/2)+L(1),x_grid(jK)-ceil(L(2)/2)+1:x_grid(jK)-ceil(L(2)/2)+L(2),3)+1;

    end
end


%moyenne des deplacement ajouter sur chaque pixel
for h1=1:Tim(1)
    for h2=1:Tim(2)
        chp_dep_int(h1,h2,1:2)=chp_dep_int(h1,h2,1:2)./chp_dep_int(h1,h2,3);
    end
end
 
chp_dep_int=chp_dep_int(:,:,1:2);


