tic

% clear all; clc;
% N = 80;
% ov = 1;                                                                                 % Oversample factor
% Ry =4; Rz =4;                                                                         % Undersampling factor
% caipi = 1;                                                                               % 1 to dephase pair lines as in 2D CAIPI
% plt = 0;                                                                                  % 1 to plot all trajectories
% sv = 0;                                                                                   % Save variables locally
% 
% rty = 1.988e-04;                                                                         % Raise time for Y grad in s
% rtz = 7.6923e-06;                                                                       % Raise time for Z grad  in s
% 
% grady = 3e-3:3e-3:30e-3;
% gradz = 40e-3:3e-3:200e-3;
% pixels = 100:50:500;
% 
% oo=0;
% info=zeros(125,5);
% 
% for nnn=1:size(pixels,2)
%     p_bw = pixels(nnn);
%     t_r = 1/p_bw;
%     limy = t_r/rty;
%     limz = t_r/rtz;
%     sinesy = 2:1:limy;
%     sinesz = 2:1:limz;
%     for jjj=1:size(grady,2)
%         gy = grady(jjj);
%         for kkk = 1:size(gradz,2)
%             gz = gradz(kkk);
%             for lll=1:size(sinesy,2)
%                 sinsy = sinesy(lll);
%                 for mmm=1:size(sinesz,2)
%                     sinsz = sinesz(mmm);

%% Add paths
%clear all; %clc
addpath(genpath('/home/amonreal/gpuNUFFT-master/gpuNUFFT'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/shepp_logan3d'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/Coilsensitivity'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/gpuNUFFT'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/general'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/Fessler_nufft'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/gpuNUFFT/gpuNUFFT-master/matlab/demo/utils'));

%% WAVE-CAIPI parameters
gy = (6e-3);                                                                            % Max amplitude of sin Y gradient
gz = (6e-3);                                                                           % Max amplitude of sin Z gradient
sinsy = 6;                                                                              % # of sins per readout line   
sinsz = 6;                                                                              % # of sins per readout line 
p_bw = 250;                                                                          % pixel BW, 70 from paper, *4 to compensate img size
gf =  1;                                                                                  % G-Factor method, 1=teoretical, 2=iterative, 3=both

%% General parameters
N = 80;
slices = 40;
ov = 6;                                                                                 % Oversample factor
Ry =4; Rz =4;                                                                         % Undersampling factor
caipi = 0;                                                                               % 1 to dephase pair lines as in 2D CAIPI
plt = 0;                                                                                  % 1 to plot all trajectories
sv = 0;                                                                                   % Save variables locally

%% Coil parameters
z_offset = [-0.025 0.025];
coil_radius = (N+30)/1000/2;                         % radius from image origin
loop_radius = N*0.06/240;                        % radius of coil
Resolution = 1e-3;

%% For Phantom
i_t = phantom3d(N);
% i_t = imnoise(i_t,'gaussian',0,0.001);                                                                
i_fov = [size(i_t,1) size(i_t,2) size(i_t,3)];                                % Vector with fov, [x y z]
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

%% Noise decorrelation
% SNR = 25;
% %add noise to data
% addpath('D:\MATLAB\External Functions\addGaussianNoise')
% %prewithen rawdata
% rawData_withNoise_Undersampled= reshape(kspace_nufft,[],nCh);
% % rawData_withNoise_Fullysampled= complex(zeros(size(kspace_nufft_fully_sampled)));
% noise_level = 0.00125* ((-6.7390) + (1i*3.1531));
% % noise_level = 0.00125*max(kspace_nufft(:));
% noise = noise_level*complex(randn(size(rawData_withNoise_Undersampled)),randn(size(rawData_withNoise_Undersampled)));
% noiseCovariance_Undersampled = noise_covariance_mtx(noise);
% rawData_withNoise_Undersampled = rawData_withNoise_Undersampled + noise;
% dmtx_Undersampled = noise_decorrelation_mtx(noiseCovariance_Undersampled);
% rawData_withNoise_Undersampled = apply_noise_decorrelation_mtx(rawData_withNoise_Undersampled,dmtx_Undersampled);
% % smaps_prew = apply_noise_decorrelation_mtx(rawData_withNoise_Undersampled,dmtx_Undersampled); 
% DECorrelatedSensitivity_Undersampled  = apply_noise_decorrelation_mtx(CoilSensitivity,dmtx_Undersampled);
% CoilSensitivity=DECorrelatedSensitivity_Undersampled;

%% Cropping image and CoilSens, if z <> x,y
i_t = i_t(:,:,((N-slices)/2)+1:((N-slices)/2)+slices);
CoilSensitivity = CoilSensitivity(:,:,((N-slices)/2)+1:((N-slices)/2)+slices,:);
CoilSensitivity = padarray(CoilSensitivity,(N*(ov-1))/2,'both');

%% Creating PSFs
[psf_y,psf_z,psf_yz]=psf(N,slices,i_fov,k_fov,ov,p_s,x_points,i_t,t,t_prime,x_calc,y_calc,z_calc,plt);

%% Calculating G-Factor
if gf == 1 || gf==3
    [g_img_t,g_av_t,g_max_t]=gfact_teo(N,slices,Ry,Rz,nCh,i_t,psf_yz,CoilSensitivity,ov);
end
if gf == 2 || gf==3
    [g_img_i,g_av_i,g_max_i]=gfact_iter(N,slices,Ry,Rz,nCh,i_t,psf_yz,CoilSensitivity,ov);
end
as(g_img_t)

%% Call next function
% Undersampling_kspace
