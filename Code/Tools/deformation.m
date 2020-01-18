function [dir deformationPrinc] = deformation(chp_dep_int,h0)

h0=h0/sum(h0);
h1=conv2([0.5;0;-0.5],1,h0);
h0=conv2([0.25;0.5;0.25],1,h0);
% displacement derivatives
dxx=conv2(h1,h0,chp_dep_int(:,:,2),'same');
dxy=conv2(h0,h1,chp_dep_int(:,:,2),'same');
dyx=conv2(h1,h0,chp_dep_int(:,:,1),'same');
dyy=conv2(h0,h1,chp_dep_int(:,:,1),'same');

for i = 1:size(chp_dep_int,1)
    for j = 1:size(chp_dep_int,2)
        
        val = [dxx(i,j) dxy(i,j); dyx(i,j) dyy(i,j)];
        % def is the Green-Lagrange deformation tensor
        def = 1/2*(val + val.' + val.'*val);

        [eigVectors eigValues]= eig (def);

        maxStrain(i,j) = max(max(eigValues));
        minStrain(i,j) = min(min(eigValues));
        if max(max(eigValues)) == 0
            dir.X(i,j)= 0;
            dir.Y(i,j)= 0;
            deformationPrinc(i,j)=max(max(eigValues));
        else
            % eigValues is a triangular matrix, hence indexX = indexY
            [indexX , indexY] = find (eigValues == max (max(eigValues)));
            % Determine principal direction of deformation i.e.
            % directions of eigen vectors associated with biggest
            % and smallest eigen values
            deformationPrinc(i,j)=max(max(eigValues));
            deformationTemp = eigVectors (:, indexX);
            dir.X(i,j)= deformationTemp(1);
            dir.Y(i,j)= deformationTemp(2);
        end

    end
end