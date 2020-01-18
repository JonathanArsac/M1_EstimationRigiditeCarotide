function writeimg(filename,data)

fid = fopen(filename,'wb');
imim = [];
for i=1:size(data,1)
    imim((i-1)*size(data,2)+1 : i*size(data,2)) = data(i,:);
end 
fwrite(fid,imim,'double');fclose(fid); 