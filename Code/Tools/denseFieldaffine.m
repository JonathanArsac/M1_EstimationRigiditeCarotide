function [chp_dep_int,pts_ax,pts_lat] = denseFieldaffine(varargin)

% dense field computation

%INPUT
% varargin{1}=chp_dep_est : estimated bilinear parameters
% varargin{2}= x_grid : lateral grid with estimated points on the initial
% image
% varargin{3}=y_grid : axial grid with estimated points on the initial
% image
% varargin{4}=Grid: mesh grid on the initial image
% varargin{5}=type_interp: 1 = moyenne; 2 = médiane

%OUTPUT
%  chp_dep_int : dense field
% pts_ax,pts_lat: mesh of the dense field


% Adrian Basarab, juin 2006


%points de calcul axial et latéral
pts_ax=varargin{3}(1)+varargin{4}(1)/2:varargin{4}(1):varargin{3}(end)-varargin{4}(1)/2;
pts_lat=varargin{2}(1)+varargin{4}(2)/2:varargin{4}(2):varargin{2}(end)-varargin{4}(2)/2;


Mpt=length(pts_ax);
Npt=length(pts_lat);
Mgd=length(varargin{3});
Ngd=length(varargin{2});

%recherche des points d'estimation les plus proches du point P
proch_ax=abs(pts_ax'*ones(1,Mgd)-ones(Mpt,1)*varargin{3});
proch_lat=abs(pts_lat'*ones(1,Ngd)-ones(Npt,1)*varargin{2});
if_ax=(proch_ax==min(proch_ax')'*ones(1,Mgd));
if_lat=(proch_lat==min(proch_lat')'*ones(1,Ngd));

for i=1:Mpt
    yy=find(if_ax(i,:));
    for k=1:numel(yy)
        gr_ax(i,k)=yy(k);
    end
end


progress_bar=timerbar(0,'calculating...');

for j=1:Npt

    timerbar((j-1)/Npt,progress_bar);

    xx=find(if_lat(j,:));
    for k=1:numel(xx)
        gr_lat(j,k)=xx(k);
    end

    for i=1:Mpt
        dep_ax=0;
        dep_lat=0;
        nb=1;
        %on considre les quatre points P (ceux ou on connait la valeur du champ estimé) autour du point considére
        for m=1:2
            for n=1:2
                %pour chacun des points P
                %calcul de la distance (x,y) au point considéré
                y=pts_ax(i)-varargin{3}(gr_ax(i,m));
                x=pts_lat(j)-varargin{2}(gr_lat(j,n));

                if((x<=varargin{4}(2)/2)&(y<=varargin{4}(1)/2))
                    %Calcul de la valeur du valeur du champ interpolé grace aux formules du déplacement bilinéaire
                    u=x*varargin{1}(gr_ax(i,m),gr_lat(j,n),1,1)+y*varargin{1}(gr_ax(i,m),gr_lat(j,n),2,1)+...
                        varargin{1}(gr_ax(i,m),gr_lat(j,n),3,1);

                    v=x*varargin{1}(gr_ax(i,m),gr_lat(j,n),1,2)+y*varargin{1}(gr_ax(i,m),gr_lat(j,n),2,2)+...
                        varargin{1}(gr_ax(i,m),gr_lat(j,n),3,2);
                    %stockage de ces valeurs dans deux tableaux
                    dep_ax(nb)=u;
                    dep_lat(nb)=v;
                    nb=nb+1;
                end
            end
        end

        switch(varargin{5})
            %calcul du champ interpolé moyen au point considéré
            case 1
                %en prenant la moyenne après avoir éliminer les valeurs qui
                %ne sont pas dans l'intervalle [mean-std mean+std]
                mu=mean(dep_ax(find(mean(dep_ax)+std(dep_ax)>=dep_ax & dep_ax>=mean(dep_ax)-std(dep_ax))));
                mv=mean(dep_lat(find(mean(dep_lat)+std(dep_lat)>=dep_lat & dep_lat>=mean(dep_lat)-std(dep_lat))));
            case 2
                %en prenant la valeur médiane
                mu=median(dep_ax);
                mv=median(dep_lat);
        end
        chp_dep_int_p(i,j,1)=mu;
        chp_dep_int_p(i,j,2)=mv;

    end
end

close(progress_bar);
%interpolation à tous les points de l'image
%grille d'interpolation
[x y]=meshgrid(pts_lat,pts_ax);
[x_new y_new]=meshgrid(pts_lat(1):pts_lat(end),pts_ax(1):pts_ax(end));
%valeurs à interpoler
u=chp_dep_int_p(:,:,1);
v=chp_dep_int_p(:,:,2);
%interpolation des composantes latérale et axiale
u_new=interp2(x,y,u,x_new,y_new);
v_new=interp2(x,y,v,x_new,y_new);
[m n]=size(u_new);
chp_dep_int=NaN*ones(m,n,2);
chp_dep_int(:,:,1)=u_new;
chp_dep_int(:,:,2)=v_new;



