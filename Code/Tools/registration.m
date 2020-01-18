function [r_corr_local,NCC_global,im1,im2_initial,im2c]=registration(im,chp_dep_int,y,x,bloc);
%recalge de l'image compressé 

%%INPUT
% im: tableau 3D avec la séquence d'images im(:,:,:) 
% chp_dep_int: champ estimé interpolé
% y,x: zone ou le recalage est calculée (y: ponits axiales, x=points latérales)
% bloc: taille du pattern pour le calcul de la corrélation


%%OUTPUT
%r_corr_local: corrélation normalisée locale entre l'image recalée et l'image initiale
%im1: image initiale
%im2c:image recalée



if ~exist('bloc','var')
    bloc = [40 20];
end


ref_y = y(1); ref_x = x(1); deb = 10;

%points de l'ancienne grille:y et x
dy=floor(y(1)-ref_y);
dx=floor(x(1)-ref_x);

%calcul de la nouvelle grille a partir de la valeur des déplacements
for i=1:length(y)
    for j=1:length(x)
        nouvo=[y(i) x(j)]+[chp_dep_int(i+dy,j+dx,1) chp_dep_int(i+dy,j+dx,2)];
        tab(i,j,:)=[nouvo(1) nouvo(2)];
    end
end

%calculs des nouvelles intensités par interpolation de l'image 2
%grille de départ
x0=x(1)-deb:x(end)+deb;
y0=y(1)-deb:y(end)+deb;
[x1 y1]=meshgrid(x0,y0);
%definition de l'image 2 aux points initials
im2i=im(floor(y0(1):y0(end)),floor(x0(1):x0(end)),2);
%grille d'arrivée
new_y=tab(:,:,1);
new_x=tab(:,:,2);
%interpolation de la 2ème image=> image recalée
im2c=interp2(x1,y1,im2i,new_x,new_y,'spline');

%on prend la meme zone sur la première image: zone qui correspond a la grille de depart
im1=im(floor(y(1):y(end)),floor(x(1):x(end)),1);
im2_initial=im(floor(y(1):y(end)),floor(x(1):x(end)),2);

%coefficient de corrélation global entre image initiale et image recalée
NCC_global=sum(sum(im1.*im2c))/sqrt(sum(sum(im1.*im1))*sum(sum(im2c.*im2c)));

% correlation locale (bloc par bloc) entre image recalée et image initiale%
it=floor(size(im1)./bloc)-1;
for i=1:it(1)
    for j=1:it(2)
        im1_bloc=im1((i-1)*bloc(1)+1:(i-1)*bloc(1)+1+bloc(1),(j-1)*bloc(2)+1:(j-1)*bloc(2)+1+bloc(2));
        im2c_bloc=im2c((i-1)*bloc(1)+1:(i-1)*bloc(1)+1+bloc(1),(j-1)*bloc(2)+1:(j-1)*bloc(2)+1+bloc(2));
        r_corr_local(i,j)=sum(sum(im1_bloc.*im2c_bloc))/sqrt(sum(sum(im1_bloc.*im1_bloc))*sum(sum(im2c_bloc.*im2c_bloc)));
    end
end

% cmd = ['Global NCC = ' num2str(NCC_global)]
% figure; imagesc(r_corr_local); colorbar; title(cmd)
