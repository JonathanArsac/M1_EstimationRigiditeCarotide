% read DICOM
dicomfile = 'IM_0031';
dicomdir = '';
[im, map] = dicomread(fullfile(dicomdir, dicomfile));
disp(size(im));
im = squeeze(im(:,:,1,:));

% write to mat
matfile = append(dicomfile, '.mat');
matdir = '';
save(fullfile(matdir, matfile), 'im', 'map');