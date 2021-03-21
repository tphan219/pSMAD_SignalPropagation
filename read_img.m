function img = read_img(fname)
% read tiff stack.
info = imfinfo(fname);
nimg = length(info);
img = imread(fname);   % first member of stack
for i = 2:nimg
    img = cat(3,img, imread(fname,i));
end


% 