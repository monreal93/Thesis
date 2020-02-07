addpath(genpath('/home/amonreal/gpuNUFFT-master/gpuNUFFT'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/shepp_logan3d'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/Coilsensitivity'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/gpuNUFFT'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/general'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/Fessler_nufft'))
addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/gpuNUFFT/gpuNUFFT-master/matlab/demo/utils'));


%% Tasks
%   - Adjust ky and kz, to make sure helix get only to 500 in both direct
%   - Check all values of Coil Sensitivity section, Resolution specially
%   - Check modif needed in NUFFT part , for non-iso img
%   - Check all paths that are needed, clean up the ones not needed
%   - Generate PSFs
%   - Fix that the sin and cos start when the readout is performed, -500Ks
%   -
%   -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
N = 60;
ov = 6;                                                                                 % Oversample factor
i_t = phantom3d(N);
i_f = fftshift(fftn(i_t));                                                                  
i_fov = [size(i_t,1) size(i_t,2) size(i_t,3)];                                % Vector with fov, [x y z]
p_s = [1e-3 1e-3 1e-3];                                                         % Vector with pixel size [x y z]
voxel_size = p_s(1)*p_s(2)*p_s(3)*.001;                                % voxel size mm3
x_points = i_fov(1)*ov;                                                             % number 6 if for oversampling, change if needed
gy_gz_amplit = 6e-3;                                                              % Max amplitude of sin readout gradients
sins = 7;                                                                                % # of sins per readout line      
p_bw = 70*4;                                                                         % pixel BW, 70 from paper, *4 to compensate img size
Ry = 3; Rz = 3;                                                                       % Undersampling factor      
k_fov = 1./p_s;                                                                        % K-space FOV
gamma = 42.58e6;                                                                % Gyromagnetic ratio
caipi = true;                                                                           % 1 to dephase pair lines as in 2D CAIPI
plt = 0;                                                                                   % 1 to plot all trajectories
nCh=32;                                                                                 % number of coils

%% Calculating parameters for Gradients
t_r = 1/p_bw;
delta_t = 1e-7;
t_prime = delta_t:delta_t:t_r;
t = 0+(t_r/x_points):t_r/x_points:t_r;
f = sins/t_r;
omg = 2*pi*f;

x_calc_points = length(t_prime);
x_calc = ((-k_fov(1)/2)+(k_fov(1)/x_calc_points)):(k_fov(1))/x_calc_points:(k_fov(1)/2);
y_calc = gamma*cumsum(gy_gz_amplit*sin(omg*t_prime)*delta_t);
z_calc = gamma*cumsum(gy_gz_amplit*cos(omg*t_prime)*delta_t);

range_y = (k_fov(1)/2)-max(y_calc);
range_z = (k_fov(1)/2)-max(z_calc);
x = interp1(t_prime,x_calc,t);
y = interp1(t_prime,y_calc,t)+range_y;
z =  interp1(t_prime,z_calc,t)+range_z;
% %%%%
figure; plot3(x,y,z)
% %%%%

% Defining spacing between helix in y and z directions
del_ky = (k_fov(2)-((max(y)-min(y))*(i_fov(2)/Ry)))/(i_fov(2)/Ry)+(max(y)-min(y));
del_kz = (k_fov(3)-((max(z)-min(z))*(i_fov(3)/Rz)))/(i_fov(3)/Rz)+(max(z)-min(z));

original_y = y; 

%%% Creating null matrices for K points
kx = zeros(i_fov(2)/Ry,i_fov(3)/Rz,x_points);
ky =kx; kz = kx;

%%  Generating K-space trajectories
for j=1:(i_fov(3)/Rz)
    for l=1:(i_fov(2)/Ry)
            kx(l,j,:) = x;
            ky(l,j,:) = y;
            kz(l,j,:) = z;
            y = y -(del_ky);
%         if (caipi == 1) &&(rem(j,2) ==0)&& (((k_fov/del_ky)-i)==1)
%            y = y -del_ky;
%             break 
%         end
    end
          z = z-(del_kz);
          y = original_y;
      if caipi
         if (rem(j,2) ==1)
                y = y - (del_ky/2);
         end
      end
end

%%  Ploting all the trajerctories
if plt == 1
    for j = 1:i_fov(2)/Ry    
        xx = squeeze(kx(j,:,:));
        yy = squeeze(ky(j,:,:));
        zz = squeeze(kz(j,:,:));
        for l=1:i_fov(3)/Rz
            plot3(xx(l,:),yy(l,:),zz(l,:));
            hold on
        end
    end
    xlabel('Kx');ylabel('Ky');zlabel('Kz');view(3)
end

%%  Interpolating
% kx = kx(:);
% ky = ky(:);
% kz = kz(:);
% 
% del_k = 1./(i_fov/1000*1);              % replace *1 with the size of the pixel...
% max_k = (1./p_s);
% 
% ksg_x = del_k(1):del_k(1):max_k(1); ksg_x = ksg_x - mean(ksg_x);
% ksg_y = del_k(2):del_k(2):max_k(2); ksg_y = ksg_y - mean(ksg_y);
% ksg_z = del_k(3):del_k(3):max_k(3); ksg_z = ksg_z - mean(ksg_z);
% 
% [ksg_X, ksg_Y, ksg_Z] = meshgrid(ksg_y,ksg_x,ksg_z);
% 
% figure; scatter3(ksg_X(:),ksg_Y(:),ksg_Z(:),5.*ones(size(ksg_X(:))),abs(i_f(:))),title('K-Space');
% xlabel('Kx');ylabel('Ky');zlabel('Kz');view(3)
% InterpolatedValues = interp3(ksg_X,ksg_Y,ksg_Z,i_f,kx,ky,kz,'cubic');
% Interpolated_kspace = reshape(InterpolatedValues,size(i_f,1)/Ry,size(i_f,2)/Rz,x_points);
% 
% u_img = ifftn(Interpolated_kspace);
% 
% as(Interpolated_kspace)
% as(u_img)
%}

%% Generate coil sensitivity
z_offset = [-0.025 0.025];
coil_radius = 0.14;
loop_radius = 0.05;
Resolution = 3e-3;
[CoilSensitivity1] = GenerateCoilSensitivity4Sim(size(i_t),[Resolution Resolution Resolution],coil_radius,loop_radius,z_offset(1),nCh./2);
[CoilSensitivity2] = GenerateCoilSensitivity4Sim(size(i_t),[Resolution Resolution Resolution],coil_radius,loop_radius,z_offset(2),nCh./2);
CoilSensitivity = cat(4,CoilSensitivity1,CoilSensitivity2);
mask = ones(size(i_t));
mask(i_t == 0) = 0;
mask = repmat(mask,1,1,1,nCh);
CoilSensitivity = CoilSensitivity.*mask;
CoilSensitivity = CoilSensitivity./1;%max(CoilSensitivity(:));

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
osf = 2; wg = 3; sw = 8;
imwidth = N;
w = 1; % no density compensation for now
FT = gpuNUFFT(tr,col(ones(size(kx(:)))),osf,wg,sw,[N,N,N],[],true);
kspace_nufft = FT*i_t;

if plt == 1
figure; view(2)
scatter3(kx(:),ky(:),kz(:),ones(size(kz(:))).*50,abs(kspace_nufft),'.')
xlabel('Kx');ylabel('Ky');zlabel('Kz');view(3)
end

kspace_nufft = reshape(kspace_nufft,20,20,360);

%% Generating PSF
% Creating cartersian trajectory
del_c = 1./(i_fov.*p_s);
x_c = zeros(size(i_t,1),size(i_t,2),x_points);
y_c = x_c;  z_c = x_c;

% xc = (k_fov(1)/2*-1)+(del_c(1)/ov):del_c(1)/ov:k_fov(1)/2;
% yc = repmat(i_fov(2)/2*-1+(p_s(2)),1,size(xc,2));
% zc = repmat(i_fov(3)/2*-1+(p_s(3)),1,size(xc,2));
xc = (k_fov(1)/2*-1)+(del_c(1)/ov):del_c(1)/ov:k_fov(1)/2;
yc = repmat((p_s(2)),1,size(xc,2));
zc = repmat((p_s(3)),1,size(xc,2));
figure; plot3(xc,yc,zc);
original_yc = yc;

for j=1:i_fov(2)
    for l=1:i_fov(3)
            x_c(l,j,:) = xc;
            y_c(l,j,:) = yc;
            z_c(l,j,:) = zc;
            yc = yc +(p_s(2));
    end
          zc = zc+(p_s(3));
          yc = original_yc;
end

%%  Ploting all the trajerctories
% if plt == 1
    for j = 1:i_fov(2)   
        for l=1:i_fov(3)
            xx = x_c(l,j,:); xx = xx(:)';
            yy = y_c(l,j,:); yy = yy(:)';
            zz = z_c(l,j,:); zz = zz(:)';
            plot3(xx,yy,zz);
            hold on
        end
    end
    xlabel('Kx');ylabel('y');zlabel('z');view(3); title('Hybrid K-space')
% end

%% Creating PSFs
psf_y = zeros(size(ky,1)*Ry,x_points);
psf_z = zeros(size(kz,1)*Rz,x_points);

t_r = 1/p_bw;
delta_t = 1e-7;
t_prime = 0:delta_t:t_r;
t = 0+(t_r/x_points):t_r/x_points:t_r;
f = sins/t_r;
omg = 2*pi*f;

x_calc_points = length(t_prime);
x_calc = ((-k_fov(1)/2)+(k_fov(1)/x_calc_points)):(k_fov(1))/x_calc_points:(k_fov(1)/2);
y_calc = gamma/2*pi*cumsum(gy_gz_amplit*sin(omg*t_prime)*delta_t);
z_calc = gamma/2*pi*cumsum(gy_gz_amplit*cos(omg*t_prime)*delta_t);
figure; plot(x_calc,y_calc)
hold on
plot(x_calc,z_calc)

Px = interp1(t_prime,x_calc,t);
Py = interp1(t_prime,y_calc,t);
Pz =  interp1(t_prime,z_calc,t);
figure; plot(Px,Py)
hold on
plot(Px,Pz)
% % %%%%
% figure; plot3(Px,Py,Pz)

% psf Y
    for l=1:i_fov(2)
        yy = y_c(l,1,:); yy = yy(:)';
        psf_y(l,:) = exp(-i*2*pi*Py.*yy);
    end

% psf Z
    for l=1:i_fov(3)
        zz = z_c(1,l,:); zz = zz(:)';
        psf_z(l,:) = exp(-i*2*pi*Pz.*zz);
    end

figure; imagesc(angle(psf_y));colormap jet, axis image off
figure; imagesc(angle(psf_z));colormap jet, axis image off

