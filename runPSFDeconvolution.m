% Example Script for PSF Deconvolution

% Load EPI testdata
load('sampledata.mat');
imageEPI = data.imageEPI;
fieldMap = data.fieldMap;

% Set parameters
params.nx = size(imageEPI,1); params.ny = size(imageEPI,2); 
params.Te = 0.034; params.Ta = 0.15; 
params.yres = 0.192*1e3/params.nx;
params.alpha = 5;
scan.type = 'epi';
t2star = 0.2;

% Reconstruct data
reconImageZeroFill = zeros(size(imageEPI));
imageConjFill = zeros(size(imageEPI));
reconImageConjFill = zeros(size(imageEPI));

w = waitbar(0,'Correcting images ...');
for Nslice = 1:size(imageEPI,3)
    waitbar(Nslice/size(imageEPI,3), w, 'Correcting images ...');
    
    % (A) Zerofilling
    scan.pf = 5/8; scan.pftype = 'zerofill'; scan.direction = 'down';
    reconImageZeroFill(:,:,Nslice) = PSF_correction_method_patzig(imageEPI(:,:,Nslice), fieldMap(:,:,Nslice), t2star*ones(size(fieldMap(:,:,Nslice))), params, scan);
    
    % (B) Conjugatefilling
    ksp_tmp = ifftshift(ifftn(fftshift(imageEPI(:,:,Nslice))));
    kspConj = rot90(conj(ksp_tmp),2);
    if strcmp(scan.direction, 'up')
        ksp_tmp(ceil(scan.pf*size(ksp_tmp,1)):end,:) = kspConj(ceil(scan.pf*size(ksp_tmp,1)):end,:);
    elseif strcmp(scan.direction, 'down')
        ksp_tmp(1:end-ceil(scan.pf*size(ksp_tmp,1)),:) = kspConj(1:end-ceil(scan.pf*size(ksp_tmp,1)),:);
    end
    imageConjFill(:,:,Nslice) = fftshift(fftn(ifftshift(ksp_tmp)));
    
    scan.pftype = 'conjfill';
    reconImageConjFill(:,:,Nslice) = PSF_correction_method_patzig(imageConjFill(:,:,Nslice), fieldMap(:,:,Nslice), t2star*ones(size(fieldMap(:,:,Nslice))), params, scan);
end
delete(w);

% Plot Results
figure
imagesc(abs([reconImageZeroFill(:,:,Nslice) reconImageConjFill(:,:,Nslice)]))
colormap(gray); axis image; axis off;

% Example for DEPICTING (not included in the testdata)
params.nx = size(imageDEPICTING,1); params.ny = size(imageDEPICTING,2);
params.Te = 0.005; params.Ta = 0.15; 
params.yres = 0.192*1e3/params.nx;
scan.type = 'depicting';
t2star = 0.2;

Nslice = 10;
reconImage{Nslice} = PSF_correction_method_patzig(imageDEPICTING(:,:,Nslice), fieldMap(:,:,Nslice), t2star*ones(size(fieldMap(:,:,Nslice))), params, scan);
