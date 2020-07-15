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

%% G-Factor method selection
gf =  5;                                                                % G-Factor method, 1=teoretical, 2=iterative, 3=fast, 4=pseudo 5=all

%% For Phantom                   
param = [];
param.N = 224;                                                    % In plane resolutio
param.slices = 120;                                               % number of slices
param.gy = (6e-3);                                              % Max amplitude of sin Y gradient
param.gz =(6e-3);                                              % Max amplitude of sin Z gradient
param.sinsy = 7;                                                 % # of sins per readout line   
param.sinsz = 7;                                                % # of sins per readout line 
param.p_bw = 70;                                             % pixel BW
param.Ry = 2;                                                     % Undersampling factor in Y 
param.Rz = 2;                                                      % Undersampling factor in Z
param.caipi_del=1;                                             % 2D-CAIPI delta
param.caipi = 0;                                                 % 1 to for 2D-CAIPI sampling
param.p_s = [1e-3 1e-3 2e-3];                            % Vector with pixel size [x y z]
param.plt = 0;                                                    % 1 to plot all trajectories
[param.ov,rg,param.if_max] = ov(param);
param.gpu = true;
param.cs = 1;                                                     % CoilSensitivity, 1=from file, 0=simulated
param.nufft = false;                                            % Use NUFFT to create Wave image, false = use psf only
i_t = phantom3d(param.N); 
param.gz_sr = 5200;                                                      % Z Gradient slew rate
param.gy_sr = 200;                                                        % Y Gradient slew rate
% Testing... REMOVEE....
param.ov = 6;

%% Calculating max sines
param = max_sines(param);

%% Other paramets
%%% Receiver coil parameters
param.nCh=32;                                                   % number of channels
param.coil_radius = (param.N+50)/1000;           % radius from image origin
param.loop_radius = param.N*0.06/240;            % radius of coil   
param.z_offset = [-3*(3/4)*param.loop_radius -(3/4)*param.loop_radius (3/4)*param.loop_radius 3*(3/4)*param.loop_radius ];
param.Resolution = 3e-3;

%%% Other parameters
param.i_fov = [param.N param.N param.slices];                                % Vector with fov, [x y z]
param.voxel_size = param.p_s(1)*param.p_s(2)*param.p_s(3);          % voxel size m3
param.x_points = (param.i_fov(1)*param.ov)-1;                                     
param.k_fov = 1./param.p_s;                                                              % K-space FOV
param.gamma = 42.58e6;                                                                 % Gyromagnetic ratio

%% Calculating parameters for Gradients
[x,y,z,x_calc,y_calc,z_calc,del_ky,del_kz,t,t_prime]=grad_param(param);

%% Cropping image and CoilSens, if z ~= x,y
if exist('info') ~= 1
    i_t = i_t(:,:,((param.N-param.slices)/2)+1:((param.N-param.slices)/2)+param.slices);
end

%% Generate coil sensitivity
CoilSensitivity = coil_sens(i_t,param);

%% Cropping image and CoilSens, if z ~= x,y
% if exist('info') ~= 1 && param.nufft
%     CoilSensitivity = CoilSensitivity(:,:,((param.N-param.slices)/2)+1:((param.N-param.slices)/2)+param.slices,:);
% end

%% Generate image w/coil sensitivity for 32 channels
i_t = repmat(i_t,1,1,1,32);
i_t = i_t.*CoilSensitivity;

%%  Generating K-space Corkscrew trajectories
[kx,ky,kz,r_ps]=cork_traj(x,y,z,del_ky,del_kz,param);
clear x y z del_ky del_kz

%% Creating PSFs
[~,~,psf_yz]=psf(i_t,t,t_prime,x_calc,y_calc,z_calc,param);
% clear x_points t t_prime x_calc y_calc z_calc

%% Wave_caipi image
%with NUFFT
if param.nufft
    % Generating NUFFT
    kspace_nufft = wc_nufft(i_t,kx,ky,kz,param);
    clear kx ky kz
    % Adding zeros to skipping positions due to 2D CAIPI and generate WAVE-CAIPI image
    i_wc = imag_wc(kspace_nufft,param);
else
    %with PSF
    i_wc = imag_wc_psf(CoilSensitivity,psf_yz,i_t,param);
end

%% Calculating G-Factor
% % Theoretical
% if gf == 1 || gf==5
%     [g_img_t,g_av_t,g_max_t]=gfact_teo1(i_t,psf_yz,CoilSensitivity,param);
% end
% % Fast
% % if gf == 3 || gf==5
% %     g_mean = gfact_iter1(psf_yz,CoilSensitivity,param);
% % end

%% SENSE Reconstruction
[i_wc_recon,rmse] = wc_sense_recon(i_wc,CoilSensitivity,psf_yz,param);
as(i_wc_recon)
