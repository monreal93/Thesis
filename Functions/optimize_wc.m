% Script to reconstruct images, it is called when using script
% "op_wc_main", it does pretty much the same as script "optimize_wc_param"

function [g_av_t,g_max_t,rmse,i_wc_recon,g_img_t,r_ps] = optimize_wc(param,ref_img)

%% Calculating max sines
param = max_sines(param);

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

%% Cropping image and CoilSens, if z ~= x,y
if param.N > param.slices
    i_t = i_t(:,:,((param.N-param.slices)/2)+1:((param.N-param.slices)/2)+param.slices);
elseif param.N < param.slices
    i_t = padarray(i_t,[0 0 (param.slices-param.N)/2],0,'both');
end

%% Generate coil sensitivity
CoilSensitivity = coil_sens(i_t,param);

%% Generate image w/coil sensitivity for 32 channels
i_t = repmat(i_t,1,1,1,32);
i_t = i_t.*CoilSensitivity;

%%  Generating K-space Corkscrew trajectories
[kx,ky,kz,r_ps]=cork_traj(x,y,z,del_ky,del_kz,param);
clear x y z del_ky del_kz

%% Creating PSFs
[~,~,psf_yz]=psf(i_t,t,t_prime,x_calc,y_calc,z_calc,param);
clear x_points t t_prime x_calc y_calc z_calc

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
% Theoretical
[g_img_t,g_av_t,g_max_t]=gfact_teo1(psf_yz,CoilSensitivity,param);

%% SENSE Reconstruction
% NUFFT undersampling
if param.nufft
    [i_wc_recon] = wc_sense_recon(i_wc,CoilSensitivity,psf_yz,param);
%     as(i_wc_recon)
else
%psf undersampling
    [i_wc_recon] = wc_sense_recon(i_wc,CoilSensitivity,psf_yz,param);
%     as(i_wc_recon)
end

%% RMSE
rmse = norm_mse(ref_img,i_wc_recon,CoilSensitivity,param);

end