% This script is to simulate a Wave-CAIPI recon for only one specific configuration
% - All the parameters can be selected in section 2), check all of them so script
%   does not stop.
% - Be careful with parameters (Nx,Ny,slices,gy,gz,sinsy,sinz) to avoid
%   that Matlab crashes with large values

%% 1) Add paths
clear all; clc
addpath(genpath('/home/amonreal/gpuNUFFT-master/gpuNUFFT'))
addpath(genpath('./Functions/Code4Alejandro/shepp_logan3d'))
addpath(genpath('./Functions/Code4Alejandro/Coilsensitivity'))
addpath(genpath('./Functions/Code4Alejandro/gpuNUFFT'))
addpath(genpath('./Functions/Code4Alejandro/general'))
addpath(genpath('./Functions/Code4Alejandro/Fessler_nufft'))
addpath(genpath('./Functions/Code4Alejandro/gpuNUFFT/gpuNUFFT-master/matlab/demo/utils'));

%% 2) Parameters                   
param = [];
param.Nx = 80;                                                    % In plane resolutio
param.Ny = 80;
param.slices = 80;                                               % number of slices
param.gy = (2e-3);                                              % Max amplitude of sin Y gradient
param.gz =(2e-3);                                              % Max amplitude of sin Z gradient
param.sinsy = 7;                                                 % # of sins per readout line   
param.sinsz = 7;                                                % # of sins per readout line 
param.p_bw = 70;                                             % pixel BW
param.Ry = 4;                                                     % Undersampling factor in Y 
param.Rz = 4;                                                      % Undersampling factor in Z
param.caipi_del=1;                                             % 2D-CAIPI delta
param.caipi = 1;                                                 % 1 to for 2D-CAIPI sampling
param.p_s = [1e-3 1e-3 1e-3];                            % Vector with pixel size [x y z]
param.plt = 0;                                                    % 1 to plot all trajectories (better to set it to 0)
[param.ov,rg,param.if_max] = ov(param);            % Calculates oversampling factor              
param.cs = 1;                                                     % CoilSensitivity, 1=from file, 0=simulated
param.nufft = false;                                            % Use NUFFT to create Wave image, false = use psf only
param.gpu = true;                                              % Use gpu, (only works if param.nufft set to true)
i_t = phantom3d(param.Nx);                               % Creating phantom
param.gz_sr = 1300;                                           % Z Gradient slew rate
param.gy_sr = 200;                                            % Y Gradient slew rate                              
ref_img = 'phantom';                                          % Reference image for RMSE, 'phantom' or 'path'
% ref_img = '/home/amonreal/Documents/Thesis/GitHub/Thesis/Data/wc_scan_150720/recon/brain_R1_1_recon.mat';

%% 3) Calculating max sines
param = max_sines(param);

%% 4) Other paramets (to simulate Coil Sensitivity and calcualted parameters)
%%% Receiver coil parameters
param.nCh=32;                                                   % number of channels
param.coil_radius = (param.Nx+50)/1000;           % radius from image origin
param.loop_radius = param.Nx*0.06/240;            % radius of coil   
param.z_offset = [-3*(3/4)*param.loop_radius -(3/4)*param.loop_radius (3/4)*param.loop_radius 3*(3/4)*param.loop_radius ];
param.Resolution = 3e-3;

%%% Calculated params
param.i_fov = [param.Nx param.Ny param.slices];                                % Vector with fov, [x y z]
param.voxel_size = param.p_s(1)*param.p_s(2)*param.p_s(3);          % voxel size m3
param.x_points = (param.i_fov(1)*param.ov)-1;                                     
param.k_fov = 1./param.p_s;                                                              % K-space FOV
param.gamma = 42.58e6;                                                                 % Gyromagnetic ratio
param.t_r = 1./param.p_bw;

%% 5) Calculating parameters for Gradients
[x,y,z,x_calc,y_calc,z_calc,del_ky,del_kz,t,t_prime]=grad_param(param);

%% 6) Cropping image and CoilSens, if z ~= x,y
if exist('info') ~= 1
    i_t = i_t(:,:,((param.Nx-param.slices)/2)+1:((param.Nx-param.slices)/2)+param.slices);
end

%% 7) Generate coil sensitivity
CoilSensitivity = coil_sens(i_t,param);

%% 8) Generate image w/coil sensitivity for 32 channels
i_t = repmat(i_t,1,1,1,32);
i_t = i_t.*CoilSensitivity;

%% 9) Generating K-space Corkscrew trajectories
[kx,ky,kz,r_ps]=cork_traj(x,y,z,del_ky,del_kz,param);
param.r_ps = r_ps;
clear x y z del_ky del_kz

%% 10) Creating PSFs
[~,~,psf_yz]=psf(i_t,t,t_prime,x_calc,y_calc,z_calc,param);
clear x_points t t_prime x_calc y_calc z_calc

%% 11) Wave_caipi image
%with NUFFT
if param.nufft
    % Generating NUFFT
    kspace_nufft = wc_nufft(i_t,kx,ky,kz,param);
%     clear kx ky kz
    % Adding zeros to skipping positions due to 2D CAIPI and generate WAVE-CAIPI image
    i_wc = imag_wc(kspace_nufft,param);
else
    %with PSF
    i_wc = imag_wc_psf(CoilSensitivity,psf_yz,i_t,param);
end

%% 12) Calculating G-Factor
% Theoretical
[g_img_t,g_av_t,g_max_t]=gfact_teo1(psf_yz,CoilSensitivity,param);
% Fast (Not sure if function works nice)
%     g_mean = gfact_iter1(psf_yz,CoilSensitivity,param);

%% 13) SENSE Reconstruction
[i_wc_recon] = wc_sense_recon(i_wc,CoilSensitivity,psf_yz,param);

%% 14) RMSE
rmse = norm_mse(ref_img,i_wc_recon,CoilSensitivity,param);

%% 15) Saving data for PSF analysis
res.i_wc_recon = i_wc_recon;
res.gf_map = g_img_t;
res.param = param;
res.rmse = rmse;
res.g_av = g_av_t;
res.g_max = g_max_t;
res.psf_yz = psf_yz;

% Saving volumes and results
% save('Data/psf_analysis/res_gd_sd.mat','res')
