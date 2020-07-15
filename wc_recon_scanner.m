%clear all; clc
clearvars -except PSFData smaps kspace phn i_wc_complete
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/PhilipsDataReaders/'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/'))
addpath(genpath('/home/jeroen/matlab/ArrShow'))
addpath(genpath('/home/jeroen/FRIDGE_jeroen/HIGHGRID/AlejandroData/Optimize_wc_4RIBSserver/functions'))

%% Parameters
param.N = 224;
param.slices = 60;
param.ov = 6;
param.sinsy = 7;
param.sinsz = 7;
param.nCh = 32;
param.p_bw = 70;
param.p_s = [1 1 2]*0.001;
param.i_fov = [param.N param.N param.slices].*param.p_s;        % [readout phase slice]
param.Ry = 1;                                                     % Undersampling factor in Y 
param.Rz = 1;                                                      % Undersampling factor in Z
param.caipi_del=1;                                             % 2D-CAIPI delta
param.caipi = 1;                                                 % 1 to for 2D-CAIPI sampling
param.x_points = (param.i_fov(1)*param.ov)-1;
param.k_fov = 1./param.p_s;
psf_source = 1;                                                 % Psf source, 1 from single projection scan (to be fitted) , 0 from field camera (loaded from file)

%% Load data
% load('./Data/wc_scan_150720/CoilSensitivitySpherePhantomNoSliceOversampling_new2mm_withMask.mat')
% load('./Data/wc_scan_150720/RawDataSphere_noSliceOversampling.mat')
% CoilSensitivity = smaps((448-224)/2:((448-224)/2)+223,:,:,:);
% CoilSensitivity = smaps;
% i_wc = kspace;

%% PSF
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
    
%     psf_z = conj(psf_z);
%     psf_y = conj(psf_y);
    
    % Masking psfs ..
    msk_y = rssq(PSFData.PSF_z_ref,4);
    msk_y(msk_y<100) = 0; msk_y(msk_y>100) = 1;
    msk_z = rssq(PSFData.PSF_y_ref,4);
    msk_z(msk_z<20) = 0; msk_z(msk_z>20) = 1;
    psf_y = msk_y.*psf_y;
    psf_z = msk_z.*psf_z;
end
    
    % cropping psfs
    if psf_source == 1
        nz_y = find(psf_y(1,:));
        nz_z = find(psf_z(1,:));
    else 
        nz_y = find(psf_y(1,:));
        nz_z = 112-60:112+59;
    end
    
    % Creating cartesian coordinates
    yc = param.i_fov(2)/2*-1:param.p_s(2):(param.i_fov(2)/2)-param.p_s(2);
    yc = yc - mean(yc);
    zc_2mm = param.i_fov(3)/2*-1:param.p_s(3):(param.i_fov(3)/2)-param.p_s(3);
    zc_2mm = zc_2mm - mean(zc_2mm);
    zc_1mm = yc;
    xc = param.i_fov(1)*param.ov/2*-1:param.p_s(1):(param.i_fov(1)*param.ov/2)-0.001;
    xc = xc - mean(xc);
    
   %shifting --> AQ
    dt = 9.300000034272795e-06;
    ds_w = 13.017600059509277e-3;
    fq = 1/(ds_w/7);
    K_ramp = 42.577e6*0.5*0.006*(0.006/200)*2*pi;
    I = round(((ds_w/7)/2)/dt);
    I = I/2;
    rmp = zc_1mm.*K_ramp;
    t_acq = 0:dt:dt*(1344-1);
    psf_z_fit = zeros(size(psf_z,1),size(psf_z,2));
    delta_t = 1e-7;
    
    psf_z_ang = angle(psf_z);
    for idx=1:size(psf_z,2)
        psf_z_tmp = psf_z_ang(:,idx);
        psf_z_tmp = unwrap(psf_z_tmp,[],1);
        
        % approach 1
%         psf_z_tmp = psf_z_tmp - mean(psf_z_tmp);
        
%         % approach 2
%         psf_z_tmp = diff(psf_z_tmp)./delta_t;
%         psf_z_tmp = circshift(psf_z_tmp,-48);
%         psf_z_tmp = cumsum(psf_z_tmp*delta_t);
%         psf_z_tmp =  interp1(psf_z_tmp,1:1344,'linear','extrap');
%         psf_z_tmp = psf_z_tmp - mean(psf_z_tmp);
% %         psf_z_tmp = psf_z_tmp.*zc_1mm(idx);               
% %         psf_z_tmp = padarray(psf_z_tmp,[1 0],'pre');

        
%         psf_z_tmp = psf_z_tmp*1.2;  %% TEstiiingggnngngn
%         psf_z_tmp = psf_z_tmp - mean(psf_z_tmp);% + min(psf_z_tmp);
%         psf_z_tmp = psf_z_tmp + 0.1763;  % TEstinggg...
        
%         psf_z_tmp = circshift(psf_z_tmp,I,1);
        
%                 % Testing...
%         amp = max(psf_y_tmp(:));
%         psf_y_sample = 42.58e6*cumsum(amp*sin(2*pi*fq*t_acq)*1e-7);
%         psf_y_sample_fit(:,idx) = exp(-1i*psf_y_sample);
        
        %    psf_y_tmp = wrapToPi(psf_y_tmp);
        %    psf_y_tmp =  exp(-1i*psf_y_tmp);
        
        if idx==75
            1;
        end
        
        psf_z_fit(:,idx) = psf_z_tmp;
    end
%     psf_z_fit = psf_z_fit - (mean(psf_z_fit)*1.2);
    
    psf_z = psf_z_fit;
    clear psf_z_fit;
    
    % Fitting psf Y and Z
    % psf_y_ang = angle(psf_y(:,nz_y));
    psf_y_ang = angle(psf_y(:,nz_y));
    % psf_z_ang = (psf_z(:,nz_z));
    
    for idx = 1:size(psf_y,1)
        yc_nz = yc(1,nz_y);
        zc_nz = zc_2mm(1,nz_z);
        
        % Unwrap
        psf_y_tmp = psf_y_ang(idx,:);
        psf_y_tmp = unwrap(psf_y_tmp,[],2);
        psf_z_tmp = psf_z(idx,:);
        psf_z_tmp = nonzeros(psf_z_tmp)';
        psf_z_tmp = unwrap(psf_z_tmp,[],2);
        
        % polyfit
        p = polyfit(yc_nz,psf_y_tmp,1);
        psf_y_tmp = polyval(p,yc);
        p = polyfit(zc_nz,psf_z_tmp,1);
        psf_z_tmp = polyval(p,zc_2mm);
        
        
        %     % Wrap back
        psf_y_tmp = wrapToPi(psf_y_tmp);
        psf_z_tmp = wrapToPi(psf_z_tmp);
        
        % Complex exponential
        psf_y_angl(idx,:) = psf_y_tmp;
        psf_z_angl(idx,:) = psf_z_tmp;
        
        psf_y_tmp =  exp(-1i*psf_y_tmp);
        psf_z_tmp =  exp(-1i*psf_z_tmp);
        
        psf_y_fit(idx,:) = psf_y_tmp;
        psf_z_fit(idx,:) = psf_z_tmp;
        
    end
    
    
    psf_y = psf_y_fit;
    psf_z = psf_z_fit;
    
    if psf_source == 1
        I = 50;  % for single projection scan
        I = 42+6;
%         I = 0;
    elseif psf_source == 0
        I = 60;  % for field camera data..
        I = 50;
    end
%     psf_z = circshift(psf_z,I,1);
    
%     psf_y = circshift(psf_y,10,1);
%[psf_y,psf_z] = deal(psf_z,psf_y);
psf_yz = repmat(psf_y,[1,1,size(psf_z,2)]) .* repmat(permute(psf_z, [1,3,2]), [1,size(psf_y,2),1]);

psf_z_unw = unwrap(angle(psf_z),[],2);
psf_y_unw = unwrap(angle(psf_y),[],2);

 %clear psf_y psf_z
 
%% Retrospectively undersampling
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
lwy = 1; upy=ceil(param.N/param.Ry);
lwz = 1; upz=ceil(param.slices/param.Rz);
i_wc = i_wc(:,lwy:upy,lwz:upz,:);

%% SENSE recon for Wave-CAIPI
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

if param.N == param.slices
    shift_amount = shift_amount .* (size(i_wc,2)/(param.Ry/param.caipi_del));
else
    shift_amount = shift_amount .* (size(i_wc,2)/param.Rz);
    shift_amount = floor(shift_amount);
end
        
for iter_wc = 1:last_iter_wc
    %tic
    sl_wc = iter_wc;
    i_wc_sl = i_wc(:,:,sl_wc,:); 
%     cs = zeros(size(i_wc_sl,1)/param.ov-1,size(i_wc_sl,2)*param.Ry-1,param.Rz,param.nCh);
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

%% Trying different things
i_wc_complete = FFT_3D_Edwin(kspace,'image');
phn = i_wc_recon;
phn(phn~=0) = 1;
phn = padarray(phn,(1344-224)/2,'both');
phn = FFT_1D_Edwin(phn,'kspace',1);
phn_wc = psf_yz.*phn;
phn_wc = FFT_1D_Edwin(phn_wc,'image',1);

% trying to reverse the wc image using the psf
wc_rev = FFT_1D_Edwin(i_wc_complete,'kspace',1);
wc_rev = wc_rev./psf_yz;
wc_rev = FFT_1D_Edwin(wc_rev,'image',1);
wc_rev = wc_rev((param.ov*param.N/2)-param.N/2:(param.ov*param.N/2)+param.N/2-1,:,:,:);
wc_rev = rssq(wc_rev,4);
as(wc_rev)

% % power spectrum of reconstructed image
% n = param.ov*param.N;
% pw_complete = FFT_1D_Edwin(i_wc_complete,'kspace',1);
% pw_complete = mean(pw_complete,4);
% pw_sample_complete = pw_complete(:,112,30);
% pw_sample_complete = abs(pw_sample_complete).^2/n;
% f = (-n/2:n/2-1)*(fq/n);
% % figure; plot(f,pw_sample_complete)
% 
% % power spectrum of reconstructed image
% n = param.ov*param.N;
% pw_complete = FFT_1D_Edwin(phn_wc,'kspace',1);
% pw_complete = mean(pw_complete,4);
% pw_sample_complete = pw_complete(:,112,30);
% pw_sample_complete = abs(pw_sample_complete).^2/n;
% f = (-n/2:n/2-1)*(fq/n);
% % figure; plot(f,pw_sample_complete)

% Making template of phn_wc and i_wc_complete to compare
phn_template = phn_wc;
phn_template(abs(phn_template)>0.1) = 1;
phn_template(abs(phn_template)<0.1) = 0;
i_wc_template = i_wc_complete;
i_wc_template(abs(i_wc_template)>3) = 5;
i_wc_template(abs(i_wc_template)<3) = 0;

% as(phn_template)
% as(i_wc_template)
as(phn_template+i_wc_template)


% psf_yz_fieldcamera = psf_yz;
% psf_yz_fieldcamera_msk = abs(FFT_1D_Edwin(psf_yz_fieldcamera,'kspace',1));
% psf_yz_fieldcamera_msk(psf_yz_fieldcamera_msk>0.1) = 1;
% psf_yz_fieldcamera_msk(psf_yz_fieldcamera_msk<0.1) = 0;
% 
% psf_yz_singleslice = psf_yz;
% psf_yz_singleslice_msk = abs(FFT_1D_Edwin(psf_yz_singleslice,'kspace',1));
% psf_yz_singleslice_msk(psf_yz_singleslice_msk>0.2) = 10;
% psf_yz_singleslice_msk(psf_yz_singleslice_msk<0.2) = 0;
% 
% as(psf_yz_fieldcamera_msk+psf_yz_singleslice_msk)


