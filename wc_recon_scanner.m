% This script reconstruct data from the scanner, it retrospectively
% under-samples the data:
% 1)    select parameters used for the image acquisition, important parameters
%       are (Ry,Rz,caipi,caipi_del) and psf_source (to check where the psf is
%       coming from)
% 2)    Load data... some lines are commented, so only run that section one
%       (it takes some time to load)
% 3)    Loads and masks PSF, if psf_source == 0 (from field camera) change
%       change the location of this data
% 4)    Pads zeros, to change size of image, so recon algorithm works.. 
% 

%% Add paths 
clear all; clc
clearvars -except PSFData smaps kspace phn i_wc_complete
addpath(genpath('./Functions/Code4Alejandro/'))

%% 1) Parameters
param.Nx = 224;
param.Ny = 224;
param.slices = 60;
param.ov = 6;
param.gy = 6e-3;
param.gz = 6e-3;
param.sinsy = 7;
param.sinsz = 7;
param.nCh = 32;
param.p_bw = 70;
param.p_s = [1 1 2]*0.001;
param.i_fov = [param.Nx param.Ny param.slices].*param.p_s;        % [readout phase slice]
param.Ry = 1;                                                     % Undersampling factor in Y 
param.Rz = 1;                                                      % Undersampling factor in Z
param.caipi_del=1;                                             % 2D-CAIPI delta
param.caipi = 0;                                                 % 1 to for 2D-CAIPI sampling
param.x_points = (param.i_fov(1)*param.ov)-1;
param.k_fov = 1./param.p_s;
psf_source = 1;                                                 % Psf source, 1 from single projection scan (to be fitted) , 0 from field camera (loaded from file)
ref_img = '/home/jeroen/FRIDGE_jeroen/HIGHGRID/AlejandroData/Data/wc_scan_150720/recon/brain_fully.mat';

%% 2) Load data
% To use brain data: replace SpherePhantom with Brain
% load('./Data/wc_scan_150720/CoilSensitivityBrainNoSliceOversampling_new2mm_withMask.mat')
% To use brain data: replace Sphere with Brain
% load('./Data/wc_scan_150720/RawDataBrain_noSliceOversampling.mat')
CoilSensitivity = smaps;
i_wc = kspace;

%% 3) Loading and masking PSF
if psf_source == 0
    % Load field camera psf, obtained from example file...
    load('Data/psf_y.mat')
    load('Data/psf_z.mat')    
elseif psf_source == 1
    %From single projection scan
%     load('./Data/wc_scan_150720/RawDataPSF.mat')
    % z  -> phase  , y -> slice
    psf_z = PSFData.PSF_y_on./PSFData.PSF_y_ref;
    psf_y = PSFData.PSF_z_on./PSFData.PSF_z_ref;
    psf_y = mean(psf_y,4);
    psf_z = mean(psf_z,4);
    
    psf_z = conj(psf_z);
    psf_y = conj(psf_y);
    
    % Masking psfs ..
    msk_y = rssq(PSFData.PSF_z_ref,4);
    msk_y(msk_y<100) = 0; msk_y(msk_y>100) = 1;
    msk_z = rssq(PSFData.PSF_y_ref,4);
    msk_z(msk_z<20) = 0; msk_z(msk_z>20) = 1;
    psf_y = msk_y.*psf_y;
    psf_z = msk_z.*psf_z;
end

%% 4) Resizing CoilSens and psf for underampling
    if param.Ry == 3
        % Padding cero to make it 240x240x60 in y direction, only or R 3x3
        CoilSensitivity = padarray(CoilSensitivity,[0 8 0 0], 'both');
%         i_wc = padarray(i_wc,[0 8 0 0], 'both');
        psf_y = padarray(psf_y,[0 8 0 0],'both');
        
    elseif param.Ry == 4
        % Padding cero to make it 224x224x64 in y direction, only or R 4x4
        CoilSensitivity = padarray(CoilSensitivity,[0 0 2 0], 'both');
        psf_z = padarray(psf_z,[0 0 2 0],'both');
    end
    
%% 5) Interpolating PSF
    % Creating cartesian coordinates
    yc = param.i_fov(2)/2*-1:param.p_s(2):(param.i_fov(2)/2)-param.p_s(2);
    yc = yc - mean(yc);
    zc_2mm = param.i_fov(3)/2*-1:param.p_s(3):(param.i_fov(3)/2)-param.p_s(3);
    zc_2mm = zc_2mm - mean(zc_2mm);
    
    for idx = 1:size(psf_y,1)
        % Finding indices of non-zero elements
        nz_y = find(msk_y(idx,:));
        nz_z = find(msk_z(idx,:));
        
        % Getting spatial coordinates of non-zero elements
        yc_nz = yc(1,nz_y);
        zc_nz = zc_2mm(1,nz_z);
       
        % Unwrap
        psf_y_tmp = angle(psf_y(idx,nz_y));
        psf_y_tmp = unwrap(psf_y_tmp,[],2);
        psf_z_tmp = angle(psf_z(idx,nz_z));
        psf_z_tmp = nonzeros(psf_z_tmp)';
        psf_z_tmp = unwrap(psf_z_tmp,[],2);
        
        % polyfit
        p = polyfit(yc_nz,psf_y_tmp,1);
        psf_y_tmp = polyval(p,yc);
        p = polyfit(zc_nz,psf_z_tmp,1);
        psf_z_tmp = polyval(p,zc_2mm);
        
        % Wrap back
        psf_y_tmp = wrapToPi(psf_y_tmp);
        psf_z_tmp = wrapToPi(psf_z_tmp);
        
        % Saving angles of psf
        psf_y_angl(idx,:) = psf_y_tmp;
        psf_z_angl(idx,:) = psf_z_tmp;
        
        % Complex exponential
        psf_y_tmp =  exp(-1i*psf_y_tmp);
        psf_z_tmp =  exp(-1i*psf_z_tmp);
        
        psf_y_fit(idx,:) = psf_y_tmp;
        psf_z_fit(idx,:) = psf_z_tmp;
    end

    psf_y = psf_y_fit;
    psf_z = psf_z_fit;

    psf_yz = repmat(psf_y,[1,1,size(psf_z,2)]) .* repmat(permute(psf_z, [1,3,2]), [1,size(psf_y,2),1]);
%     clear psf_y psf_z
   
%% 6) Retrospectively undersampling
% Resizing wave-caipi image for underampling
    if param.Ry == 3
        % Padding cero to make it 240x240x60 in y direction, only or R 3x3
        i_wc = FFT_3D_Edwin(i_wc,'image');
        i_wc = padarray(i_wc,[0 8 0 0], 'both');
        i_wc = FFT_3D_Edwin(i_wc,'kspace');
    elseif param.Ry == 4
        % Padding cero to make it 224x224x64 in y direction, only or R 4x4
        i_wc = FFT_3D_Edwin(i_wc,'image');
        i_wc = padarray(i_wc,[0 0 2 0], 'both');
        i_wc = FFT_3D_Edwin(i_wc,'kspace');
    end

% Creating mask for CAIPI undersampling
msk = zeros(size(i_wc));
kz_i=1;
ky_i=1;
        for ii=1:param.Rz
                msk(:,ky_i:param.Ry*param.Rz:end,kz_i:param.Rz:end,:) = 1;
                kz_i = kz_i+param.caipi_del;
                ky_i = ky_i+param.Ry;
            if kz_i>param.Rz
                kz_i=1;
            end
        end
 % Undersampling in CAIPI way
i_wc = i_wc.*msk;
i_wc = FFT_3D_Edwin(i_wc,'image');

% Cropping from start
lwy = 1; upy=floor(param.Ny/param.Ry);
lwz = 1; upz=floor(param.slices/param.Rz);
i_wc = i_wc(:,lwy:upy,lwz:upz,:);

%% 7) SENSE recon for Wave-CAIPI
y_skip = size(i_wc, 2);
i_wc_recon = complex(zeros(size(CoilSensitivity(:,:,:,1))));
param.psf_len = size(i_wc, 1);
param.img_len = size(i_wc,1)/param.ov;
param.pad_size = (param.psf_len - param.img_len) / 2;

% Reconstructing all phantom or only one slice of i_wc
last_iter_wc = size(i_wc,3);

% Getting shift amount
shift_amount = zeros(1,param.Rz);
for ii=1:param.Rz-1
    shift_amount(ii+1) = ii;
end

if param.Nx == param.slices
    shift_amount = shift_amount .* (size(i_wc,2)/(param.Ry/param.caipi_del));
else
    shift_amount = shift_amount .* (size(i_wc,2)/param.Rz);
    shift_amount = ceil(shift_amount);
end
        
for iter_wc = 1:last_iter_wc
    %tic
    sl_wc = iter_wc;
    i_wc_sl = i_wc(:,:,sl_wc,:); 
    cs = zeros(size(CoilSensitivity,1),size(CoilSensitivity,2),param.Rz,param.nCh);
    slice_ind = zeros(1,param.Rz);
    sl_i_t = sl_wc;

    % Extracting CS of only overlaped slices and getting slice indices
    cs(:,:,1,:) = CoilSensitivity(:,:,sl_i_t,:);
    slice_ind(1) = sl_i_t;

    for ii=2:param.Rz
            if (sl_i_t+size(i_wc,3)) > size(i_wc,3)*param.Rz
                sl_i_t = sl_i_t+size(i_wc,3)-size(CoilSensitivity,3);
                cs(:,:,ii,:) = CoilSensitivity(:,:,sl_i_t,:);
                slice_ind(ii)=sl_i_t;
            else
                sl_i_t = sl_i_t+size(i_wc,3);
                cs(:,:,ii,:) = CoilSensitivity(:,:,sl_i_t,:);
                slice_ind(ii)=sl_i_t;
            end
    end
    
        msk_roi = (cs~=0);
        msk_roi = msk_roi(:,:,:,1);
        Img_WAVE = zeros(size(i_wc_sl,1)/(param.ov),size(i_wc_sl,2)*(param.Ry),size(i_wc_sl,3)*(param.Rz));
        lsqr_iter = 200;
        lsqr_tol = 1e-5; 
        
        psf_use = repmat(psf_yz(:,:,slice_ind), [1,1,1,param.nCh]);   
        receive_use = zeros(size(cs));
        mask_use = zeros(size(msk_roi));
        for nz = 1:param.Rz
            receive_use(:,:,nz,:) = circshift(cs(:,:,nz,:), [0,shift_amount(nz),0,0]);
            psf_use(:,:,nz,:) = circshift(psf_use(:,:,nz,:), [0,shift_amount(nz),0,0]);
            mask_use(:,:,nz) = circshift(msk_roi(:,:,nz), [0,shift_amount(nz),0]);
        end
        mask_use = sum(mask_use,3);

        for cey = 1:y_skip
            cey_ind = cey : y_skip : cey + (param.Ry-1) * y_skip;
            if sum(sum(mask_use(:,cey_ind),1),2) > 0
                psfs = psf_use(:,cey_ind,:,:);
                rcv = receive_use(:,cey_ind,:,:);
                rhs = squeeze( i_wc_sl(:,cey,:,:) );
                param.psfs = psfs;
                param.rcv = rcv;
                [res, flag, relres, iter] = lsqr(@apply_wave, rhs(:), lsqr_tol, lsqr_iter, [], [], [], param);        
                Img_WAVE(:,cey_ind,:) = reshape(res, [param.img_len, param.Ry, param.Rz]);            
            end
        end

        for nz = 1:param.Rz
            Img_WAVE(:,:,nz) = circshift(Img_WAVE(:,:,nz), [0,-shift_amount(nz),0,0]);
        end
        
        i_wc_recon(:,:,slice_ind) = Img_WAVE;
        
        %toc
end

%% 8) Calculating G-Factor
% Theoretical
[g_img_t,g_av_t,g_max_t]=gfact_teo1(psf_yz,CoilSensitivity,param);

%% 9) RMSE
rmse = norm_mse(ref_img,i_wc_recon,CoilSensitivity,param);

res = [];
res.i_wc_recon = i_wc_recon;
res.gf_map = g_img_t;
res.param = param;
res.rmse = rmse;
res.g_av = g_av_t;
res.g_max = g_max_t;

%% 10) Save data
% tit = sprintf('Data/wc_scan_150720/recon/%i_%i_%i_%i_%i_%i_%i_%i_%i.mat',...
%     param.Ry,param.Rz,param.Nx,param.slices,param.gy.*1000,param.gz*1000,...
%     param.sinsy,param.sinsz,param.p_bw);
% save(tit,'res')