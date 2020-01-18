function H=DeformMesh(im,patchSize,x,y,dep_ax,dep_lat,x_grid,y_grid,L)

% see example in lancer_seq.m in BDBM directory
% im = image sequence
% patchSize


N = size(im,3);

x = x(1):patchSize(1):x(end); y = y(1):patchSize(2):y(end);

for k = 1:(N-1)

    imagesc(im(:,:,k)); colormap gray;

    for B1=1:(size(x,2)-1) % pour chaque bloc affiché
        for B2=1:(size(y,2)-1)
            d1(k,B1,B2)=mean2(dep_lat((y(B2)-L(2)/2-y_grid(1)):(y(B2)+L(2)/2-y_grid(1)),(x(B1)-L(1)/2-x_grid(1)):(x(B1)+L(1)/2-x_grid(1)),k));
            d2(k,B1,B2)=mean2(dep_ax((y(B2)-L(2)/2-y_grid(1)):(y(B2)+L(2)/2-y_grid(1)),(x(B1)-L(1)/2-x_grid(1)):(x(B1)+L(1)/2-x_grid(1)),k));
        end
    end
    for B1=2:(size(x,2)-2)
        for B2=2:(size(y,2)-2)
            line([x(B1)+d1(k,B1,B2) x(B1)+patchSize(1)+d1(k,B1+1,B2)], [y(B2)+d2(k,B1,B2) y(B2)+d1(k,B1+1,B2)],'Color','w','LineWidth',0.5)          
            line([x(B1)+d1(k,B1,B2) x(B1)+d1(k,B1,B2-1)], [y(B2)+d2(k,B1,B2) y(B2)-patchSize(2)+d1(k,B1,B2-1)],'Color','w','LineWidth',0.5)          
            line([x(B1)+d1(k,B1,B2) x(B1)+d1(k,B1,B2+1)], [y(B2)+d2(k,B1,B2) y(B2)+patchSize(2)+d1(k,B1,B2+1)],'Color','w','LineWidth',0.5)
        end
    end
    H(k)=getframe
end

