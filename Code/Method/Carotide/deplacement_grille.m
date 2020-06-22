in = "/home/foxiler/git/M1_EstimationRigiditeCarotide/Donnees/save_estimation/estimation_IM_0020.mat";
f = open(in);
estimation = f.estimation;
% clearvars f

chp_dep_est = estimation.chp_dep_est;
x_grid = estimation.x_grid;
y_grid = estimation.y_grid;
x_dep_grid = estimation.x_dep_grid;
y_dep_grid = estimation.y_dep_grid;
deplacement_x = estimation.deplacement_x;
deplacement_y = estimation.deplacement_y;
frames = estimation.frames;
% clearvars estimation

load IM_0020.mat
im = double(im);

Grid=[10 10];
type_interp = 1;

pas = 5;

nb_frames = length(frames);
x_min = x_dep_grid(1);
x_max = x_dep_grid(end);
x_size = x_max - x_min + 1;
y_min = y_dep_grid(1);
y_max = y_dep_grid(end);
y_size = y_max - y_min + 1;

x_ind = 1:pas:x_size;
y_ind = 1:pas:y_size;

[X_base, Y_base] = meshgrid(x_min:x_max, y_min:y_max);
[X, Y] = meshgrid(x_min:pas:x_max, y_min:pas:y_max);



for k = 1:nb_frames
    imagesc([x_min, x_max], [y_min, y_max], im(y_min:y_max,x_min:x_max,frames(k))); colormap(gray);
    dep_x = interp2(X_base, Y_base, deplacement_x(:, :, k), X, Y);
    dep_y = interp2(X_base, Y_base, deplacement_y(:, :, k), X, Y);
    X = X + dep_x;
    Y = Y + dep_y;
    l1 = line(X, Y);
    l2 = line(transpose(X), transpose(Y));
    pause(0.1)
end

% x = 360;
% y = 190;
% 
% for k = 1:nb_frames
%     imagesc(im(:,:,frames(k))); colormap(gray);
%     dep_x = interp2(X_base, Y_base, deplacement_(:, :, k), x, y);
%     dep_y = interp2(X_base, Y_base, deplacement_y(:, :, k), x, y);
%     x = x + dep_x;
%     y = y + dep_y;
%     hold on
%     scatter(x,y);
%     hold off
%     pause(0.1);
% end