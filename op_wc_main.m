% This script will simulate different Wave-CAIPI configurations, it will call the function "optimize_wc" to simulate the reconstruction.
% Section 2) is to selct parameters for the recon
% Section 3) is to select the WAVE-CAIPI parameters to be tested
% Section 6) performs the simulations


%% 1) Add paths
clear all; clc;
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/PhilipsDataReaders/'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/'))
addpath(genpath('/home/jeroen/matlab/ArrShow'))
addpath(genpath('/home/jeroen/FRIDGE_jeroen/HIGHGRID/AlejandroData/Optimize_wc_4RIBSserver/functions'))

%% 2) Select prameters and simulation options
Nx=80;
Ny=80;
slices=80;
Ry=4;
Rz=4;
caipi_del=2;                                                     % CAIPI Delta, Note: not sure if RECON works for all values
gz_sr = 5200;                                                   % Z Gradient slew rate
gy_sr = 200;                                                     % Y Gradient slew rate
gpu = false;                                                      % Use gpu for NUFFT
cs = 1;                                                              % CoilSensitivity, 1=from file, 0=simulated
nufft = false;                                                     % Use NUFFT to create Wave image, false = use psf only
ref_img = 'phantom';                                        % Reference image for RMSE, 'phantom' or 'path'
% ref_img = '/home/amonreal/Documents/Thesis/GitHub/Thesis/Data/wc_scan_150720/recon/brain_R1_1_recon.mat';
sv_data = false;                            % Save all volume reconstructions???
warning('off','all')

%% 3) Manually setting parameters to try
gy = [0.010 0.020 0.030];
gz = [0.010 0.030 0.050];
sinsy = [6 10 20];
sinsz = [10 20 40];
p_bw = [100,200,400];

%% 4) Loop over different values
i_data = 1;
for ip_bw=1:length(p_bw)
    for igy=1:length(gy)
        for igz=1:length(gz)
                    % Calculating max sines
                    t_r = 1./p_bw(ip_bw);
                    gz_rt = 1/(gz_sr/gz(igz));           % Z Gradient raise time
                    gy_rt = 1/(gy_sr/gy(igy));          % Y Gradient raise time

                    %%
            for isinsy=1:length(sinsy)
                for isinsz=1:length(sinsz)
                    tic
                    %% 5) Parameters into structure
                    param = [];
                    param.gf =  1;                                                  % G-Factor method, 1=teoretical, 2=iterative, 3=fast, 4=pseudo 5=all
                    param.Nx = Nx;                                                    % In plane resolution
                    param.Ny = Ny;                                                    % In plane resolution
                    param.slices =slices;                                         % number of slices
                    param.gy = gy(igy);                                           % Max amplitude of sin Y gradient
                    param.gz = gz(igz);                                            % Max amplitude of sin Z gradient
                    param.sinsy = sinsy(isinsy);                               % # of sins per readout line
                    param.sinsz = sinsz(isinsz);                                % # of sins per readout line
                    param.p_bw = p_bw(ip_bw);                               % pixel BW
                    param.Ry = Ry;                                                   % Undersampling factor in Y
                    param.Rz = Rz;                                                    % Undersampling factor in Z
                    param.caipi_del=caipi_del;                                             % 2D-CAIPI delta
                    param.caipi = 1;                                                 % 1 to for 2D-CAIPI sampling
                    param.p_s = [1e-3 1e-3 1e-3];                            % Vector with pixel size [x y z]
                    param.plt = 0;                                                     % 1 to plot all trajectories
                    param.cs = cs;
                    param.nufft = nufft;
                    param.gpu = gpu;
                    param.gz_sr = gz_sr;
                    param.gy_sr = gy_sr;
                    param.t_r = 1/param.p_bw;
                    ov_max = floor((t_r/(1e-6))./N);
                    [param.ov,rg,param.if_max] = ov(param);
                    if param.ov > ov_max
                        param.ov = ov_max;
                    end
                    
                    %% 6) Simulate Wave-CAIPI
                    [g_av,g_max,rmse,i_wc_recon,gf_map,r_ps] = optimize_wc(param,ref_img);
                    
                    %% 7) Saving data
                    data(i_data,1) = param.gy;
                    data(i_data,2) = param.gz;
                    data(i_data,3) = param.sinsy;
                    data(i_data,4) = param.sinsz;
                    data(i_data,5) = param.p_bw;
                    data(i_data,6) = g_av;
                    data(i_data,7) = g_max;
                    data(i_data,8) = rmse;
                    data(i_data,9) = r_ps(1);
                    data(i_data,10) = r_ps(2);
                    data(i_data,11) = r_ps(3);
                    data(i_data,12) = param.ov;
                    i_data = i_data+1;
                    fprintf('Done with %d out of %d.\n',i_data-1,size(gy,2)*size(gz,2)*size(sinsy,2)*size(sinsz,2)*size(p_bw,2));
                    
                    % Saving reconstruction volumes
                    if sv_data
                        res = [];
                        res.i_wc_recon = i_wc_recon;
                        res.gf_map = gf_map;
                        res.param = param;
                        res.rmse = rmse;
                        res.g_av = g_av;
                        res.g_max = g_max;
                        tit = sprintf('Data/R_%i_%i_recon/%i_%i_%i_%i_%i_%i_%i_%i_%i.mat',...
                            param.Ry,param.Rz,param.Nx,param.slices,param.gy.*1000,param.gz*1000,...
                            param.sinsy,param.sinsz,param.p_bw,param.Ry,param.Rz);
                        save(tit,'res')
                    end
                    toc
                end
            end
        end
    end
end

%% 8) Saving data file..
file = sprintf('Data/Data_%i_%i_%i_%i_%i_%i_%i_%i_%i.mat',...
    param.Nx,param.slices,param.gy*1000,param.gz*1000,param.sinsy,param.sinsz,param.p_bw,param.Ry,param.Rz);
save(file,'data');
