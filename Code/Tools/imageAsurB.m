function [B0,A0,mask] = imageAsurB(A,B,x,y,MinMax,val_t,type_affichage)
%affichage de l'image A en premier plan sur l'image B en arriere plan
%[B0,A0,mask] = imageAsurB(A,B,x,y,[MinMax],[val_t],[type_affichage]);
%INPUT
%  A :   image en premier plan
%  B :   image en arriere plan
%  x :   abscisses de la zone d'insertion de B dans A
%  y :   ordonnees de la zone d'insertion de B dans A
%  MinMax :   dynamique de l'image A fixee a [Min Max]
%  val_t :   valeur transparence entre 0 et 1 (1 : opaque)
%  type_affichage :   'bord_parchemin' ou 'bord_net' ou 'giboulee'
%OUTPUT
%  B0 :   image en arriere plan recadree
%  A0 :   image en premier plan recadree
%  mask : masque de sur-impression

%  Ph. Delachartre, decembre 2004

% valeurs par defaut
if ~exist('MinMax','var')
    MinMax = []; 
end
if ~exist('val_t','var')
    val_t = 0.7; 
end
if ~exist('type_affichage','var')
    type_affichage = 'bord_parchemin'; 
end

%représentation des valeurs de l'image A entre +/- K sigma
%moy = mean2(A);
%sigma = std2(A);
%Aout = find(abs(A- moy) > K*sigma);
%if ~isempty(Aout), A(Aout) = 0; end

% représentation des valeurs de l'image dans l'intervalle [Min Max]
if ~isempty(MinMax), 
    Aout = find((A < MinMax(1)) | (A > MinMax(2))); 
    if ~isempty(Aout), A(Aout) = 0; end
end

% recadrage de l'image B par rapport a A
B0 = recadrage(B,A);

% rendu du bord de l'image
maskA = renduBordImage(A,x,y,val_t,type_affichage);

%%% affichage de A sur B %%%
mask = ones(size(B));
mask(y,x) = maskA;

A0 = zeros(size(B));
A0(y,x) = A;
i = find(A0==0);
if ~isempty(i), mask(i) = 0; end

figure
imagesc(B0);
hold on
h = imagesc(A0);
set(h,'AlphaData',mask);
ylim([1 size(B0,1)])
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          FONCTIONS INTERNES        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B0 = recadrage(B,A)
% recadrage de l'image B par rapport a A
A_min = min(A(:));
A_max = max(A(:));
B_min = min(B(:)); 
B_max = max(B(:)); 
a = (A_max-A_min)/(B_max-B_min);
b = A_max-a*B_min;
B0 = a*B + b;
return

function mask = renduBordImage(A,x,y,val_t,type_affichage)
% rendu du bord de l'image
Ly = length(y);
Lx = length(x);
Lh = floor(0.1*size(A))+1;
Lh = Lh + rem(Lh,2);
Lh2 = Lh/2;

switch type_affichage
   case 'bord_parchemin'
       mask1 = (rand(size(A)) > 0.5)&(rand(size(A)) > 0.5)...
       &(rand(size(A)) > 0.5);
       %filtre gaussien
       nh1 = [1:Lh(1)]-Lh2(1);sh1 = floor(Lh(1)/5)+1;
       nh2 = [1:Lh(2)]-Lh2(2);sh2 = floor(Lh(2)/8)+1;
       h = exp(-(nh1'/sh1).^2)*exp(-(nh2/sh2).^2);
       
       mask2 = filter2(h/sum(h(:)),mask1);
       mask = (mask2 > mean2(mask2))*val_t;
       mask(1+Lh2(1):Ly-Lh2(1), 1+Lh2(2):Lx-Lh2(2)) = ones(size(A)-Lh)*val_t;
	case 'bord_net'
        mask = val_t*ones(size(A)); 
        
    case 'giboulee'
       mask1 = (rand(size(A)) > 0.5)&(rand(size(A)) > 0.5)...
       &(rand(size(A)) > 0.5);
       se = strel('arbitrary',[0 1 0;1 1 1;1 1 1;0 1 0]);
       MN = floor(Lh./[8 8])+1;
       %se = strel('rectangle',MN);
       se = strel('arbitrary',(rand(MN) > 0.5));
       mask = imdilate(mask1,se);
       %mask2 = mask1;
       %mask2(1+Lh2(1):Ly-Lh2(1), 1+Lh2(2):Lx-Lh2(2)) = ones(size(S)-Lh);

   otherwise
      error('option inconnue')
end

return

