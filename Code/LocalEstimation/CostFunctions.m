function [best_corr_pos, val_c]=CostFunctions(varargin);

% function to estimate the 2D spatial delay

% %varargin{1}=PA : block to search
% varargin{2}=ZR : search zone
% varargin{3}=vI : interpolation factor
% varargin{4}= type_correl :type of local 2D spatial delay estimation.
%                                           1 = SAD ; 2 = normalized correlation; 3 = SAD + gradient descend search;
%                                           4 = elimination successive;


a=varargin{1};
L=size(a);                  %taille pattern
iR=size(varargin{2})-(L-1).*varargin{3};       %it�ration sur la zone de recherche
if varargin{4}==2
    indice= isnan(varargin{1});       %correction pixel manquant
    varargin{1}(find(indice==1))=0;
end
val_c=inf;

switch(varargin{4})
    case 1 %recherche exhastive
        for i=1:iR(1)
            for j=1:iR(2)
                PAP=varargin{2}(i:varargin{3}(1):i+(L(1)-1)*varargin{3}(1),j:varargin{3}(2):j+(L(2)-1)*varargin{3}(2));  %pattern test
                r_corr= sum_diff (varargin{1},PAP);
                tab_c(i,j)=r_corr;
                if r_corr<val_c                         % si la valeur est inf�rieure
                    best_corr_pos=[i j];                % nouveau point de r�f�rence
                    val_c=r_corr;
                end
            end
        end

    case 2 %corr�lation normalis� - non interpol�
        for i=1:iR(1)
            for j=1:iR(2)
                PAP=varargin{2}(i:varargin{3}(1):i+(L(1)-1)*varargin{3}(1),j:varargin{3}(2):j+(L(2)-1)*varargin{3}(2));  %pattern test
                r_corr=1-(sum(sum(varargin{1}.*PAP))/sqrt(sum(sum(varargin{1}.*varargin{1}))*sum(sum(PAP.*PAP)))); %corr�lation normalis�e
                tab_c(i,j)=r_corr;
                if r_corr<val_c                         % si la valeur est inf�rieure
                    best_corr_pos=[i j];                % nouveau point de r�f�rence
                    val_c=r_corr;
                end
            end
        end

    case 3 %descente de gradient -non interpol�

        %tableau de marqueur permmetant de savoir si le pixel consid�r� a �t� visit�: si le marqueur est �gale � 1, le pixel n'a jamais �t� encore visit�.
        marqueur=ones(iR(1),iR(2));
        %position initial pour la descente de gradient
        [pos_min,min]=init(iR,L,varargin{3},varargin{2},varargin{1});
        nb=0;
        lim_inf=[1 1]; %sert a delimiter la zone de recherche
        while(nb~=8) %on continue tant que l'on trouve un voisin ayant une corr�lation inf�rieure � la pr�c�dente
            nb=0;
            %on se postione sur le dernier pixel trouv� ayant la plus petite corr�lation
            %on cherche si un de ses voisins pr�sente une plus petite corr�lation
            pos_centre=pos_min;
            for j=-1:1
                for i=-1:1
                    %on prend un des huit voisins
                    pos=[pos_centre(1)+j pos_centre(2)+i];
                    %si ce pixel existe, c'est-�-dire si il appartient � la zone de recherche (pour les cas ou on se trouve sur les bords de la zone de recherche)
                    if(pos>=lim_inf & pos<=iR)
                        %si ce pixel n'est pas le centre et qu'il n'a pas encore �t� visit�
                        if(( pos(1)~=pos_centre(1)|pos(2)~=pos_centre(2))& marqueur(pos(1),pos(2))==1)
                            %calcul du pattern correpondant
                            PAP=varargin{2}(pos(1):varargin{3}(1):pos(1)+(L(1)-1)*varargin{3}(1),pos(2):varargin{3}(2):pos(2)+(L(2)-1)*varargin{3}(2));
                            %calcul du MAD entre ce nouveau pattern et le pattern de r�f�rence
                            MAD=sum_diff(varargin{1},PAP);
                            %comparaison  entre MAD et le minimum de corr�lation
                            if(MAD<min) %si MAD est plus petit, le pixel consid�r� devient le nouveau minimum
                                marqueur(pos_min(1),pos_min(2))=0;
                                %on change la valeur du minimum
                                pos_min=pos;
                                min=MAD;
                            else        %sinon le minimum reste inchang�
                                marqueur(pos(1),pos(2))=0;
                                nb=nb+1;    %incr�mentation de la variable nb car ce voisin est non valide (= il a �t� visit� et sa corr�lation n'est pas inf�rieur au minimum)
                            end
                        elseif(( pos(1)~=pos_centre(1)|pos(2)~=pos_centre(2))& marqueur(pos(1),pos(2))==0)  %si le voisin a deja �t� visit�e et qu'il ne soit pas le centre
                            nb=nb+1;
                        end
                    else %si le pixel n'appartient pas a la zone de recherche, il est consid�r� comme un pixel non valide
                        nb=nb+1;
                    end
                end
            end
        end %on s'arrete d�s que nb=8, c'est-�-dire quand tous les voisins sont non valides
        best_corr_pos=pos_min;
        val_c=min;

    
    case 4 %m�thode de l'�limination successive
        %pour les cas ou les deux images sont interpol�s
        %interpolation du pattern de r�f�rence
        [x y]=meshgrid(1:L(2),1:L(1));
        [x_new y_new]=meshgrid(1:(1/varargin{3}(2)):L(2),1:(1/varargin{3}(1)):L(1));
        PA_inter=interp2(x,y,varargin{1},x_new,y_new);
        %taille du pattern et nombre d'it�rations correspondante
        L_inter=size(PA_inter);
        iR_inter=size(varargin{2})-(L_inter-1);
        %somme des valeurs abolus de la pattern de ref�rence interpol� PA_inter
        R=sum(sum(abs(PA_inter)));
        %calcul de la somme de la norme pour les patterns interpol�s de la zone de recherche
        SN=sum_norm(varargin{2},L_inter);
        %initialisation du premier pattern interpol�
        m=1;n=1;
        PAP_inter2=varargin{2}(m:m+L_inter(1)-1,n:n+L_inter(2)-1);
        %calcul de MAD initial
        MAD=sum_diff(PA_inter,PAP_inter2);
        %calcul de la somme de la norme pour ce pattern initial
        M=SN(m,n);
        for i=1:iR_inter(1)
            for j=1:iR_inter(2)
                %construction du pattern interpol� � la position i,j
                PAP_inter=varargin{2}(i:i+L_inter(1)-1,j:j+L_inter(2)-1);
                %calcul de la somme de la norme pour ce pattern via le tableau SN
                M=SN(i,j);
                %si la valeur satisfait � cette relation , on considere qu'il est assez semblabe au pattern de r�f�rence donc on le prend en compte
                if(R-M<=MAD & M-R<=MAD)
                    MAD1=sum_diff(PA_inter,PAP_inter);
                    if(MAD1<MAD) %si MAD1 est inf�rieur � MAD,cette position devient la nouvelle r�f�rence
                        MAD=MAD1;
                        best_corr_pos=[i j];
                        val_c=MAD;
                    end
                end
            end
        end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fonctions anexes pour calculer les fonctions de co�t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [MAD]= sum_diff (PA,PAP)
%calcul SAD: somme de la valeur absolue de la diff�rence entre le pattern de r�f�rence PA et le pattern candidat PAP
Pdiff=PA-PAP;
indice= isnan(Pdiff);       %correction pixel manquant
Pdiff(find(indice==1))=0;
MAD=sum(sum(abs(Pdiff)));;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [pos_min,val_min]=init(iR,L,vI,ZR,PA)

% fonction pour determiner la position de d�part pour la descente de gradient

%points initiaux possibles
pos=zeros(5,2);
pos(1,:)=2;
pos(2,:)=[2 iR(2)-1];
pos(3,:)=[floor(iR(1)/2) floor(iR(2)/2)];
pos(4,:)=[iR(1)-1 2];
pos(5,:)=iR(2)-1;

%calcul des patterns et de la corr�lation pour chaque point pour chaque point
for m=1:5
    PAP1=ZR(pos(m,1):vI(1):pos(m,1)+(L(1)-1)*vI(1),pos(m,2):vI(2):pos(m,2)+(L(2)-1)*vI(2));
    val(m)=sum_diff(PA,PAP1);
end

[val_min pos_val]=min(val);
pos_min=pos(pos_val,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function[SN]=sum_norm(ZR,L)

%initialisation
[H W]=size(ZR);
N=L(1);
M=L(2);
R1=zeros(N,W);
C=zeros(H-N+1,W);
SN=zeros(H-N+1,W-M+1);
ZR=abs(ZR);

%premiere bande de ligne
R1=ZR(1:1+(N-1),1:W);
%% etape 1 %%%%%
% on decoupe l'image en (H-N) bandes de lignes
%on calcule la somme des colonnes pour les autres bandes
for i=1:(H-N+1)
    if(i==1)
        C(i,:)=sum(R1);
    else
        C(i,:)=C(i-1,:)-ZR(i-1,:)+ZR(i-1+N,:);
    end
end

%%%% etape 2 %%%%%
%on fait la somme pour chaque bloc de colonne de chaque bloc de lignes
% => on fait la somme des normes de chaque pattern de taille N*M de la zone de recherche
for i=1:(H-N+1)
    for j=1:(W-M+1)
        if (j==1)
            SN(i,j)=sum(C(i,j:j-1+M));
        else
            SN(i,j)=SN(i,j-1)-C(i,j-1)+C(i,j-1+M);
        end
    end
end
