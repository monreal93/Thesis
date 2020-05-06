%% Tasks
%   - Check modif needed in NUFFT part , for non-iso img
%   - 
%   - 
%   - 
%   - 
%   -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add paths
clear all; clc
addpath(genpath('/home/amonreal/gpuNUFFT-master/gpuNUFFT'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/shepp_logan3d'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/Coilsensitivity'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/gpuNUFFT'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/general'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/Fessler_nufft'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/gpuNUFFT/gpuNUFFT-master/matlab/demo/utils'));
addpath(genpath('/home/amonreal/Documents/Thesis/GitHub/Thesis/Functions'));

%%  For NIFTI
% %%% User input:
% file = 'maria.nii.gz';
% % file = 'sub-OAS30001_ses-d0757_run-01_T1w.nii.gz';
% i_t = niftiread(file); info = niftiinfo(file);
% i_f = fftshift(fftn(i_t));
% i_fov = [size(i_t,1) size(i_t,2) size(i_t,3)];                                   % mm3
% p_s = [info.PixelDimensions(1) info.PixelDimensions(2) info.PixelDimensions(3)]*0.001;
% voxel_size = p_s(1)*p_s(2)*p_s(3)*.001;
% x_points = i_fov(1)*6;
% gy_gz_amplit = 6e-3;
% sins = 7;
% p_bw = 70;
% Ry = 4; Rz = 2;
% k_fov = 1./p_s;
% gamma = 42.58e6;                                                                % Gyromagnetic ratio
% caipi = 1;
% % 1 to dephase pair lines as in 2D CAIPIA
% plt = 0;                                                                                   % 1 to plot all trajectories

%% For Phantom
% %%% User input:                   
N = 80;                                                                               % In plane resolution
slices =80;                                                                          % number of slices
ov =4;                                                                                 % Oversample factor
gy = (6e-3);                                                          % Max amplitude of sin Y gradient
gz = (6e-3);                                                          % Max amplitude of sin Z gradient
sinsy = 6;                                                                              % # of sins per readout line   
sinsz = 6;                                                                              % # of sins per readout line 
p_bw = 450;                                                                          % pixel BW, 70 from paper, *4 to compensate img size
Ry = 2; Rz =2;                                                                         % Undersampling factor
caipi = 1;                                                                               % 1 to dephase pair lines as in 2D CAIPI
caipi_del=1;
plt = 0;                                                                                  % 1 to plot all trajectories
sv = 0;                                                                                   % Save variables locally
% %%%

%%% Receiver coil parameters
z_offset = [-0.025 0.025];
coil_radius = (N+30)/1000/2;                         % radius from image origin
loop_radius = N*0.06/240;                        % radius of coil
Resolution = 1e-3;

%%% Other parameters
i_t = phantom3d(N); 
% msk=i_t;
% msk(msk~=0) = 1;
% i_t = imnoise(i_t,'gaussian',0,0.001);
% i_t = i_t.*msk;
i_fov = [N N slices];                                % Vector with fov, [x y z]
p_s = [1e-3 1e-3 1e-3];                                                         % Vector with pixel size [x y z]
voxel_size = p_s(1)*p_s(2)*p_s(3);                                         % voxel size m3
x_points = (i_fov(1)*ov)-1;                                                      % number 6 if for oversampling, change if needed    
k_fov = 1./p_s;                                                                       % K-space FOV
gamma = 42.58e6;                                                               % Gyromagnetic ratio
nCh=32;                                                                                % number of coils

%% Calculating parameters for Gradients
[x,y,z,x_calc,y_calc,z_calc,del_ky,del_kz,t,t_prime]=grad_param(gy,gz,sinsy,sinsz,p_bw,p_s,k_fov,i_fov,Ry,Rz,x_points,gamma);

%% Generate coil sensitivity
CoilSensitivity = coil_sens(z_offset,coil_radius,loop_radius,Resolution,i_t,nCh,plt);

%% Cropping image and CoilSens, if z <> x,y
i_t = i_t(:,:,((N-slices)/2)+1:((N-slices)/2)+slices);
CoilSensitivity = CoilSensitivity(:,:,((N-slices)/2)+1:((N-slices)/2)+slices,:);

%% Generate K-space w/coil sensitivity for 32 channels
i_t = repmat(i_t,1,1,1,32);
i_t = i_t.*CoilSensitivity;

%%  Generating K-space Corkscrew trajectories
[kx,ky,kz]=cork_traj(x,y,z,del_ky,del_kz,i_fov,Ry,Rz,x_points,caipi,plt);
clear x y z del_ky del_kz

%% Generating NUFFT
kspace_nufft = wc_nufft(kx,ky,kz,k_fov,N,slices,i_t,ov,Ry,Rz,nCh,plt);
clear kx ky kz

%% Noise decorrelation
% SNR = 25;
% %add noise to data
% addpath('D:\MATLAB\External Functions\addGaussianNoise')
% %prewithen rawdata
% rawData_withNoise_Undersampled= reshape(kspace_nufft,[],nCh);
% % rawData_withNoise_Fullysampled= complex(zeros(size(kspace_nufft_fully_sampled)));
% noise_level = 0.00125*max(kspace_nufft(:));
% noise = noise_level*complex(randn(size(rawData_withNoise_Undersampled)),randn(size(rawData_withNoise_Undersampled)));
% noiseCovariance_Undersampled = noise_covariance_mtx(noise);
% rawData_withNoise_Undersampled = rawData_withNoise_Undersampled + noise;
% dmtx_Undersampled = noise_decorrelation_mtx(noiseCovariance_Undersampled);
% rawData_withNoise_Undersampled = apply_noise_decorrelation_mtx(rawData_withNoise_Undersampled,dmtx_Undersampled);
% % smaps_prew = apply_noise_decorrelation_mtx(rawData_withNoise_Undersampled,dmtx_Undersampled); 
% DECorrelatedSensitivity_Undersampled  = apply_noise_decorrelation_mtx(CoilSensitivity,dmtx_Undersampled);
% CoilSensitivity=DECorrelatedSensitivity_Undersampled;

%% Adding zeros to skipping positions due to 2D CAIPI and generate WAVE-CAIPI image
i_wc = imag_wc(N,slices,Ry,Rz,caipi,caipi_del,kspace_nufft,plt);
clear kspace_nufft 
as(i_wc)

%% Creating PSFs
[psf_y,psf_z,psf_yz]=psf(N,slices,i_fov,k_fov,ov,p_s,x_points,i_t,t,t_prime,x_calc,y_calc,z_calc,plt);
clear x_points t t_prime x_calc y_calc z_calc

%% Saving WAVE-CAPI, PSF and Sensitivity map and parameters
param = [];
param.psf_len = size(i_wc, 1);                                                      % psf length
param.img_len = size(i_wc,1)/ov;                                                % image length (without oversampling)
param.pad_size = (param.psf_len - param.img_len) / 2;               % zero pad size (psf length - image length) / 2
param.num_chan = size(i_wc,4);
param.Ry = Ry;
param.Rz = Rz;
param.ov = ov;

if sv
    save('Data/param.mat','param')
    save('Data/sens_map.mat','CoilSensitivity')
    save('Data/img_wc.mat','i_wc')
    save('Data/psf_yz.mat','psf_yz')
    save('Data/phantom.mat','i_t')
end

wc_sense_recon