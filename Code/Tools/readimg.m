function im = readimg(filename,dimx,dimy)

fid = fopen(filename,'rb');
img = fread(fid,'double');
for i=1:dimx
   im(i,:) = img((i-1)*dimy+1 : i*dimy);
end;
fclose(fid);
