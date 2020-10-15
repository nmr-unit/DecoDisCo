%% Example Script for PSF Deconvolution

%% Example 1: Partial Fourier 5/8 in-vivo
% Load (partial Fourier) EPI testdata
load('EPI_partialFourier0p56_inVivo.mat');
kspData  = data.kspData;
fieldMap = data.fieldMap;

% Set parameters
params.nx    = size(kspData,1); 
params.ny    = size(kspData,2); 
params.Te    = 0.034; % in seconds
params.Ta    = 0.20; 
params.yres  = 192/params.nx;
params.alpha = 2;
scan.type    = 'ge-epi';
t2star       = 0.05;

% Reconstruct data
imageDistorted     = zeros(size(kspData));
reconImageZeroFill = zeros(size(kspData));
imageConjFill      = zeros(size(kspData));
reconImageConjFill = zeros(size(kspData));

w = waitbar(0,'Correcting images ...');
for Nslice = 1:size(kspData,3)
    waitbar(Nslice/size(kspData,3), w, 'Correcting images ...');
    
    %% (A) Zerofilling
    % Transform k-space data to image domain
    imageDistorted(:,:,Nslice) = fftshift(fftn(ifftshift(kspData(:,:,Nslice))));
    % Set parameters
    scan.pf         = 0.56; 
    scan.pftype     = 'zerofill'; 
    scan.direction  = 'down';
    % Correct image
    reconImageZeroFill(:,:,Nslice) = PSF_correction_method_patzig(imageDistorted(:,:,Nslice), fieldMap(:,:,Nslice), t2star*ones(size(fieldMap(:,:,Nslice))), params, scan);
    
    %% (B) Conjugatefilling
    ksp_tmp = kspData(:,:,Nslice);
    % Conjugate k-space
    kspConj = rot90(conj(kspData(:,:,Nslice)),2);
    % Fill non-acquired k-space lines according to scan direction
    if strcmp(scan.direction, 'up')
        ksp_tmp(ceil(scan.pf*size(ksp_tmp,1)):end,:)   = kspConj(ceil(scan.pf*size(ksp_tmp,1)):end,:);
    elseif strcmp(scan.direction, 'down')
        ksp_tmp(1:end-ceil(scan.pf*size(ksp_tmp,1)),:) = kspConj(1:end-ceil(scan.pf*size(ksp_tmp,1)),:);
    end
    % Transform k-space data to image domain
    imageConjFill(:,:,Nslice) = fftshift(fftn(ifftshift(ksp_tmp)));
    % Set parameters
    scan.pftype = 'conjfill';
    % Correct image
    reconImageConjFill(:,:,Nslice) = PSF_correction_method_patzig(imageConjFill(:,:,Nslice), fieldMap(:,:,Nslice), t2star*ones(size(fieldMap(:,:,Nslice))), params, scan);
end
delete(w);

% Plot Results
figure('Name', ['EPI correction of slice ' num2str(Nslice)])
subplot(1,3,1); imagesc(abs(imageDistorted(:,:,Nslice))); title('Distorted image'); axis image; axis off;
subplot(1,3,2); imagesc(abs(reconImageZeroFill(:,:,Nslice))); title('Corrected image (zero filling)'); axis image; axis off;
subplot(1,3,3); imagesc(abs(reconImageConjFill(:,:,Nslice))); title('Corrected image (conjugate filling)'); axis image; axis off;
colormap(gray);

%% Example 2: Partial Fourier 5/8 in-vivo
% Load (partial Fourier) EPI testdata
load('EPI_doubleCenterOut_DEPICTING_phantom.mat');
kspData   = data.kspData;
fieldMap  = data.fieldMap;
t2starMap = data.t2starMap;

% Set parameters
params.nx    = size(kspData,1); 
params.ny    = size(kspData,2); 
params.Te    = 0.0041; % in seconds
params.Ta    = 0.096; 
params.yres  = 210/params.nx;
params.alpha = 0.05;
scan.type    = 'depicting';

% Reconstruct data
imageDistorted              = fftshift(fftn(ifftshift(kspData)));
imageCorrectedWithoutT2star = PSF_correction_method_patzig(imageDistorted, fieldMap, ones(size(t2starMap)), params, scan);
imageCorrected              = PSF_correction_method_patzig(imageDistorted, fieldMap, t2starMap, params, scan);

% Plot Results
figure('Name', 'DEPICTING correction')
subplot(1,3,1); imagesc(abs(imageDistorted)); title('Distorted image'); axis image; axis off;
subplot(1,3,2); imagesc(abs(imageCorrectedWithoutT2star)); title('Corrected image without T2*'); axis image; axis off;
subplot(1,3,3); imagesc(abs(imageCorrected)); title('Corrected image including T2*'); axis image; axis off;
colormap(gray);
