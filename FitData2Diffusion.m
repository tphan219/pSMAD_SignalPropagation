function FitData2Diffusion
%
% If make source binary, ON/OFF the support of source around one clump of secreting cells
% will set overall intensity in that area, so really a separately adjustable parameter.

data_dir = '/Users/siggia/Dropbox/Desktop/phan/image_proc/spread_test_data/';
dapi = imread([data_dir,'Timelapse_41_w1DAPI.TIF']);
smad = imread([data_dir,'Timelapse_41_w3Cy5.TIF']);
bmp  = imread([data_dir,'Timelapse_41_w2Rhodamine-RFP.TIF']);

% otsu thresholding for nuclei.
nucl = im2bw(dapi, graythresh(dapi));

% zero diffusing material outside of colony.
colony = imclose(nucl, strel('disk',10,0));
colony = imfill(colony,'holes');

% some of bmp image is saturated at 4095. make it binary, and get rid of fragments
source = bmp > 2500;
source = imopen(source, strel('disk',4,0));

% exactly zero non nuclear smad1 to cleanly distinguish background from real signal
smad(~nucl) = 0;

phi = zeros(size(source));
phi1 = diffuse_phi(phi, colony, source, 0.002, 1000);

% pull out one active region, since each region has its own scale.  Fits about the same
% with lam <~ 0.001 and niter >~ 400
phi_sub = phi1(256:512, 1:384);
smad_sub = smad(256:512, 1:384);

% get rid of non nuclei points, need also to cutoff phi for high values since smad is
% bounded (and also potentially non linear fn of BMP)
phi_plot = smad_sub > 100;  
figure, plot(phi_sub(phi_plot), smad_sub(phi_plot), '.');
return