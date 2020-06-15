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
% param = [];
% file = 'maria_ss.nii.gz';
% param.ov =3;                                                      % Oversample factor
% param.gy = (6e-3);                                              % Max amplitude of sin Y gradient
% param.gz = (6e-3);                                              % Max amplitude of sin Z gradient
% param.sinsy = 6;                                                 % # of sins in Y per readout line   
% param.sinsz = 6;                                                 % # of sins in Z per readout line 
% param.p_bw = 450;                                             % pixel BW
% param.Ry = 4;                                                     % Undersampling factor in Y
% param.Rz =4;                                                      % Undersampling factor in Z
% param.caipi_del=2;                                             % 2D-CAIPI delta
% param.caipi = 1;                                                 % 1 to for 2D-CAIPI sampling
% param.plt = 0;                                                    % 1 to plot all trajectories
% i_t = niftiread(file);
% info = niftiinfo(file);
% param.i_fov = [size(i_t,1) size(i_t,2) size(i_t,3)];
% param.p_s = [info.PixelDimensions(1) info.PixelDimensions(2) info.PixelDimensions(3)]*0.001;
% [param.N,N_ind]=max(param.i_fov);
% [param.slices,sli_ind]=min(param.i_fov);
% i_t = permute(i_t,[N_ind N_ind+1 sli_ind]);
% clear sli_ind N_ind file 

%% For Phantom                   
param = [];
param.N = 80;                                                    % In plane resolution
param.slices =80;                                               % number of slices
param.gy = (6e-3);                                              % Max amplitude of sin Y gradient
param.gz = (40e-3);                                              % Max amplitude of sin Z gradient
param.sinsy = 6;                                                 % # of sins per readout line   
param.sinsz = 40;                                                 % # of sins per readout line 
param.p_bw = 250;                                             % pixel BW
param.ov = round(((max(param.gy,param.gz)/0.03)/param.N)*sqrt((param.N^2)+(param.slices^2))+2);
param.ov = 9;
param.Ry = 4;                                                     % Undersampling factor in Y 
param.Rz = 4;                                                      % Undersampling factor in Z
param.caipi_del=2;                                             % 2D-CAIPI delta
param.caipi = 1;                                                 % 1 to for 2D-CAIPI sampling
param.p_s = [1e-3 1e-3 1e-3];                            % Vector with pixel size [x y z]
param.plt = 0;                                                    % 1 to plot all trajectories
i_t = phantom3d(param.N); 

%% Other paramets
%%% Receiver coil parameters
param.nCh=32;                                                   % number of channels
param.coil_radius = (param.N+50)/1000;           % radius from image origin
param.loop_radius = param.N*0.06/240;            % radius of coil   
param.z_offset = [-2*(3/4)*param.loop_radius 0 2*(3/4)*param.loop_radius 2*2*(3/4)*param.loop_radius ];
param.Resolution = 3e-3;

%%% Other parameters
param.i_fov = [param.N param.N param.slices];                                % Vector with fov, [x y z]
param.voxel_size = param.p_s(1)*param.p_s(2)*param.p_s(3);          % voxel size m3
param.x_points = (param.i_fov(1)*param.ov)-1;                                     
param.k_fov = 1./param.p_s;                                                              % K-space FOV
param.gamma = 42.58e6;                                                                 % Gyromagnetic ratio

%% Calculating parameters for Gradients
[x,y,z,x_calc,y_calc,z_calc,del_ky,del_kz,t,t_prime]=grad_param(param);

%% Generate coil sensitivity
CoilSensitivity = coil_sens(i_t,param);

%% Cropping image and CoilSens, if z ~= x,y
if exist('info') == 0
    i_t = i_t(:,:,((param.N-param.slices)/2)+1:((param.N-param.slices)/2)+param.slices);
    CoilSensitivity = CoilSensitivity(:,:,((param.N-param.slices)/2)+1:((param.N-param.slices)/2)+param.slices,:);
end

%% Generate image w/coil sensitivity for 32 channels
i_t = repmat(i_t,1,1,1,32);
i_t = i_t.*CoilSensitivity;

%%  Generating K-space Corkscrew trajectories
[kx,ky,kz]=cork_traj(x,y,z,del_ky,del_kz,param);
clear x y z del_ky del_kz

%% Generating NUFFT
kspace_nufft = wc_nufft(i_t,kx,ky,kz,param);
clear kx ky kz

%% Adding zeros to skipping positions due to 2D CAIPI and generate WAVE-CAIPI image
i_wc = imag_wc(kspace_nufft,param);
clear kspace_nufft
as(i_wc)

%% Creating PSFs
[~,~,psf_yz]=psf(i_t,t,t_prime,x_calc,y_calc,z_calc,param);
clear x_points t t_prime x_calc y_calc z_calc

%% SENSE Reconstruction
[i_wc_recon,rmse] = wc_sense_recon(i_wc,CoilSensitivity,psf_yz,param);
as(i_wc_recon)
