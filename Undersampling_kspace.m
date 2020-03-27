tic
%% Tasks
%   - Check modif needed in NUFFT part , for non-iso img
%   - 
%   - 
%   - 
%   - 
%   -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add paths
clear all; %clc
addpath(genpath('/home/amonreal/gpuNUFFT-master/gpuNUFFT'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/shepp_logan3d'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/Coilsensitivity'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/gpuNUFFT'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/general'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/Fessler_nufft'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/gpuNUFFT/gpuNUFFT-master/matlab/demo/utils'));

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
%%% User input:
N = 80;
ov = 6;                                                                                 % Oversample factor
gy_gz_amplit = (12e-3);                                                          % Max amplitude of sin readout gradients 
sins = 12;                                                                              % # of sins per readout line      
p_bw = 500;                                                                          % pixel BW, 70 from paper, *4 to compensate img size
Ry =4; Rz =4;                                                                         % Undersampling factor
caipi = 1;                                                                               % 1 to dephase pair lines as in 2D CAIPI
plt = 0;                                                                                  % 1 to plot all trajectories
sv = 0;                                                                                   % Save variables locally
%%%

i_t = phantom3d(N);
% i_t = imnoise(i_t,'gaussian',0,0.001);
i_f = fftshift(fftn(i_t));                                                                  
i_fov = [size(i_t,1) size(i_t,2) size(i_t,3)];                                % Vector with fov, [x y z]
p_s = [1e-3 1e-3 1e-3];                                                         % Vector with pixel size [x y z]
voxel_size = p_s(1)*p_s(2)*p_s(3);                                         % voxel size m3
x_points = (i_fov(1)*ov)-1;                                                      % number 6 if for oversampling, change if needed    
k_fov = 1./p_s;                                                                       % K-space FOV
gamma = 42.58e6;                                                               % Gyromagnetic ratio
nCh=32;                                                                                % number of coils

%% Calculating parameters for Gradients
t_r = 1/p_bw;
delta_t = 1e-7;
t_prime = 0:delta_t:t_r;
t = 0:t_r/x_points:t_r;
f = sins/t_r;
omg = 2*pi*f;
x_calc_points = length(t_prime);

%For helix....
x_calc = ((-k_fov(1)/2)+(k_fov(1)/x_calc_points)):(k_fov(1))/x_calc_points:(k_fov(1)/2);
x_calc = x_calc-mean(x_calc);
y_calc = gamma*cumsum(gy_gz_amplit*cos(omg*t_prime)*delta_t);
%y_calc = y_calc - mean(y_calc);
z_calc = gamma*cumsum(gy_gz_amplit*sin(omg*t_prime)*delta_t);
%z_calc = z_calc - mean(z_calc);
range_y = (k_fov(1)/2)-max(y_calc);
range_z = (k_fov(1)/2)-max(z_calc);
x = interp1(t_prime,x_calc,t,'linear','extrap');
y = interp1(t_prime,y_calc,t,'linear','extrap')+range_y;
z =  interp1(t_prime,z_calc,t,'linear','extrap')+range_z;

% Defining spacing between helix in y and z directions
del_ky = (1/(i_fov(2)*p_s(2))*Ry);
del_kz = (1/(i_fov(3)*p_s(3))*Rz);
original_y = y;
original_z = z;

% Creating null matrices for K points
kx = zeros(x_points+1,i_fov(2)/Ry,i_fov(3)/Rz);
ky =kx; kz = kx;

%% Generate coil sensitivity
z_offset = [-0.025 0.025];
coil_radius = 0.14;                         % radius from image origin
loop_radius = 0.08;                        % radius of coil
Resolution = 3e-3;
[CoilSensitivity1] = GenerateCoilSensitivity4Sim(size(i_t),[Resolution Resolution Resolution],coil_radius,loop_radius,z_offset(1),nCh./2);
[CoilSensitivity2] = GenerateCoilSensitivity4Sim(size(i_t),[Resolution Resolution Resolution],coil_radius,loop_radius,z_offset(2),nCh./2);
CoilSensitivity = cat(4,CoilSensitivity1,CoilSensitivity2);
mask = ones(size(i_t));
mask(i_t == 0) = 0;
mask = repmat(mask,1,1,1,nCh);
CoilSensitivity = CoilSensitivity.*mask;
CoilSensitivity = CoilSensitivity./max(CoilSensitivity(:));
if plt == 1
    as(CoilSensitivity)
end

%% Generate K-space w/coil sensitivity for 32 channels
i_t = repmat(i_t,1,1,1,32);
i_t = i_t.*CoilSensitivity;

%%  Generating K-space Corkscrew trajectories
%%% Moving in y direction first and then z...
% for j=1:(i_fov(3)/Rz)
%     for l=1:(i_fov(2)/Ry)
%             kx(:,l,j) = x;
%             ky(:,l,j) = y;
%             kz(:,l,j) = z;
%             y = y -(del_ky);
% %         if (caipi == 1) &&(rem(j,2) ==0)&& (((k_fov/del_ky)-i)==1)
% %            y = y -del_ky;
% %             break 
% %         end
%     end
%           z = z-(del_kz);
%           y = original_y;
%       if caipi
%          if (rem(j,2) ==1)
%                 y = y - (del_ky/2);
%          end
%       end
% end

%%% Different approach moving in z direction first and then y....
y = y-(del_ky/2);
z = z-(del_kz/2);
original_z = z;
for j=1:(i_fov(2)/Ry)
    for l=1:(i_fov(3)/Rz)
            kx(:,j,l) = x;
            ky(:,j,l) = y;
            kz(:,j,l) = z;
            z = z -(del_kz);
%         if (caipi == 1) &&(rem(j,2) ==0)&& (((k_fov/del_ky)-i)==1)
%            y = y -del_ky;
%             break 
%         end
    end
          y = y-(del_ky);
          z = original_z;
      if caipi
         if (rem(j,2) ==1)
                z = z - (del_kz/2);
         end
      end
end

% figure;scatter(ky(:),kz(:)); xlabel('Ky');ylabel('Kz');axis image

%  Ploting all the trajerctories
if plt == 1
 figure;
    for j = 1:i_fov(2)/Rz   
        xx = squeeze(kx(:,:,j));
        yy = squeeze(ky(:,:,j));
        zz = squeeze(kz(:,:,j));
        for l=1:i_fov(3)/Ry
            plot3(xx(:,l),yy(:,l),zz(:,l));
            hold on
        end
    end
    xlabel('Kx');ylabel('Ky');zlabel('Kz');view([90 0 0])
end

%% Generating NUFFT
% Normalize trajectories
tr = zeros(3,length(kx(:)));
tr(1,:) = kx(:)./k_fov(1);
tr(2,:) = ky(:)./k_fov(2);
tr(3,:) = kz(:)./k_fov(3);

disp('Generate NUFFT Operator without coil sensitivities');
disp_slice=N/2;
useGPU = true;
useMultiCoil = 0;
osf = 2; wg = 3; sw = 5;
imwidth = N;

i_t = padarray(i_t,(((size(i_t,2)*ov)-size(i_t,2))/2),'both');
FT = gpuNUFFT(tr,col(ones(size(kx(:)))),osf,wg,sw,[N*ov,N,N],[],true);
kspace_nufft = FT*i_t;
test_nufft = kspace_nufft;

if plt == 1
    figure; view(2)
    scatter3(kx(:),ky(:),kz(:),ones(size(kz(:))).*50,abs(kspace_nufft(:,1)),'.')
    xlabel('Kx');ylabel('Ky');zlabel('Kz'); view([90 0 0])
end
kspace_nufft = reshape(kspace_nufft,N*ov,N/Ry,N/Rz,nCh);
kspace_nufft = flip(kspace_nufft,2);
kspace_nufft = flip(kspace_nufft,3);

%% Adding zeros to skipping positions due to 2D CAIPI and generate WAVE-CAIPI image
if caipi && Ry ~= 1 && Rz ~= 1    
    kspace_new = zeros([size(kspace_nufft,1) size(kspace_nufft,2)*Ry size(kspace_nufft,3)*Rz size(kspace_nufft,4)]);
    kk=1;
 % Define the start position for the shifted row in CAIPI 
    if Ry==4
      mm=3;
    else
      mm=2;
    end
 
    for iter = 1:N/Ry
            if mod(iter,2) == 0
                kspace_new(:,kk,mm:Ry:end,:) = kspace_nufft(:,iter,:,:);
                kk = kk+Ry;
            else
                 kspace_new(:,kk,1:Ry:end,:) = kspace_nufft(:,iter,:,:);
                 kk = kk+Ry;
            end 
    end
    
    lwy = fix(((size(kspace_new,2)-size(kspace_nufft,2))/2)+1);
    upy = lwy +(size(kspace_nufft,2))-1;
    lwz = fix(((size(kspace_new,3)-size(kspace_nufft,3))/2)+1);
    upz = lwz+(size(kspace_nufft,3))-1;  

    test_img = FFT_3D_Edwin(kspace_new,'image' );
    test_img1 = test_img(:,lwy:upy,lwz:upz,:);
    i_wc = test_img1;
else
    i_wc = FFT_3D_Edwin(kspace_nufft,'image' );
end

%% Creating cartersian trajectory
del_c = 1./(i_fov.*p_s);
x_c = zeros(x_points+1,size(i_t,2),size(i_t,3));
y_c = x_c;  z_c = x_c;

xc = (k_fov(1)/2*-1)+((del_c(1)/ov)/2):del_c(1)/ov:k_fov(1)/2;
yc = repmat((i_fov(2)/2*-1),1,size(xc,2))*p_s(2)+(p_s(2)/2);
zc = repmat((i_fov(3)/2*-1),1,size(xc,2))*p_s(3)+(p_s(3)/2);
original_yc = yc;

for j=1:i_fov(3)
    for l=1:i_fov(2)
            x_c(:,l,j) = xc;
            y_c(:,l,j) = yc;
            z_c(:,l,j) = zc;
            yc = yc +(p_s(2));
    end
          zc = zc+(p_s(3));
          yc = original_yc;
end

% figure;scatter(y_c(:),z_c(:)); xlabel('Ky');ylabel('Kz');axis image

% Ploting all the trajerctories
if plt == 1
figure;
    for j = 1:i_fov(2)   
        for l=1:i_fov(3)
            xx = x_c(:,l,j); xx = xx(:)';
            yy = y_c(:,l,j); yy = yy(:)';
            zz = z_c(:,l,j); zz = zz(:)';
            plot3(xx,yy,zz);
            hold on
        end
    end
    xlabel('Kx');ylabel('y');zlabel('z');view(3); title('Hybrid K-space');
end

%% Creating PSFs
psf_y = zeros(x_points+1,size(ky,2)*Ry);
psf_z = zeros(x_points+1,size(kz,3)*Rz);

Px = interp1(t_prime,x_calc,t,'linear','extrap');
Py = interp1(t_prime,y_calc,t,'linear','extrap');
Pz =  interp1(t_prime,z_calc,t,'linear','extrap');

if plt==1
    figure; plot(Px,Py)
    hold on
    plot(Px,Pz)
end

% psf Y
    for l=1:i_fov(2)
        yy = y_c(:,l,1); yy = yy(:)';
        psf_y(:,l) = exp(-1i*2*pi*Py.*yy);
    end
    
% psf Z
    for l=1:i_fov(3)
        zz = z_c(:,1,l); zz = zz(:)';
        psf_z(:,l) = exp(-1i*2*pi*Pz.*zz);
    end 

 psf_yz = repmat(psf_y,[1,1,size(psf_z,2)]) .* repmat(permute(psf_z, [1,3,2]), [1,size(psf_y,2),1]);
 
if plt == 1
    figure; imagesc(angle(psf_y).');colormap jet, axis image off ; title('PSF Y')
    figure; imagesc(angle(psf_z).');colormap jet, axis image off; title('PSF Z')
end

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

toc