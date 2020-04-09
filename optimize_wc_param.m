tic

% clear all; clc;
% grad = 3e-3:3e-3:15e-3;
% sines = 4:4:20;
% pixels = 100:100:500;
% oo=0;
% info=zeros(125,5);
% 
% for ll=1:size(grad,2)
%     for jj=1:size(sines,2)
%             for kk=1:size(pixels,2)
%                 
% gy_gz_amplit = grad(ll);                                                          % Max amplitude of sin readout gradients 
% sins = sines(jj);                                                                      % # of sins per readout line      
% p_bw = pixels(kk);                                                                  % pixel BW, 70 from paper, *4 to compensate img size

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

%% WAVE-CAIPI parameters
gy = (6e-3);                                                          % Max amplitude of sin Y gradient
gz = (6e-3);                                                          % Max amplitude of sin Z gradient
sinsy = 6;                                                                              % # of sins per readout line   
sinsz = 6;                                                                             % # of sins per readout line 
p_bw = 800;                                                                          % pixel BW, 70 from paper, *4 to compensate img size

%% For Phantom
%%% User input:
N = 80;
ov = 1;                                                                                 % Oversample factor
Ry =2; Rz =2;                                                                         % Undersampling factor
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
fy = sinsy/t_r;
fz = sinsz/t_r;
omgy = 2*pi*fy;
omgz = 2*pi*fz;
x_calc_points = length(t_prime);

%For helix....
x_calc = ((-k_fov(1)/2)+(k_fov(1)/x_calc_points)):(k_fov(1))/x_calc_points:(k_fov(1)/2);
x_calc = x_calc-mean(x_calc);
y_calc = gamma*cumsum(gy*cos(omgy*t_prime)*delta_t);
%y_calc = y_calc - mean(y_calc);
z_calc = gamma*cumsum(gz*sin(omgz*t_prime)*delta_t);
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

% load 'Data/CoilSens_decorr_80.mat'

IntensityCorrection = 1./sqrt(sum(abs(CoilSensitivity).^2,4));
coilSensNew = CoilSensitivity.*repmat(IntensityCorrection,1,1,1,nCh);
coilSensNew(isinf(coilSensNew)) = 0;
coilSensNew(isnan(coilSensNew)) = 0;

CoilSensitivity = coilSensNew;

% CoilSensitivity=coilSensNew;
% load 'Data/CoilSens_decorr_80.mat'
% CoilSensitivity=imnoise(CoilSensitivity,'Gaussian',0,0.001);

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

%% Calculating g-factor Theoretical Method
tic

psf = repmat(psf_yz,[1 1 1 nCh]);
% aa = (FFT_1D_Edwin(coilSensNew,'kspace',1));
aa = (FFT_1D_Edwin(CoilSensitivity,'kspace',1));
enc = psf.*aa;
enc = FFT_1D_Edwin(enc,'image',1);
% enc = CoilSensitivity;

% Only collapsed colums and slices...
g_img_t = zeros(N,N,N);
for jj=1:N/Rz  %2D
sli = jj:N/Rz:N;
    for ii=1:N/Ry
        col=ii:N/Ry:N;
        enc_res = single(reshape(permute(enc(:,col,sli,: ),[1 4 2 3]),N*nCh,Ry*Rz));
        for kk=1:N
            enc_res1 = enc(kk,col,sli,:);
            enc_res1 = reshape(permute(enc_res1,[1 4 2 3]),nCh,[]);
            EE = enc_res1'*enc_res1;
            % Approach 1 Theoretical
            EE_inv = pinv((EE));
            EE_diag = diag(EE);
            EE_inv_diag = abs(diag(EE_inv));
            g_f = sqrt(EE_diag.*EE_inv_diag);
            g_f = reshape(g_f,1,Ry,Rz);
            g_img_t(kk,col,sli)=g_f;
%             % Approach 2 iterative
%             g_f = zeros(1,size(EE,1));
%             for ll=1:size(EE,1)
%                 e_rho = zeros(size(EE,1),1);
%                 e_rho(ll) = 1;
%                 [d,flag]=lsqr(EE,e_rho);
%                 g_f(ll) = sqrt(e_rho'*d)*sqrt(e_rho'*EE*e_rho);
%             end
%                 g_f = reshape(g_f,1,Ry,Rz);
%                 g_img(kk,col,sli)=g_f;
        end
    end
end

disp('Theoretical method: ')
toc

msk=i_t;
msk(msk~=0) = 1;
g_img_t = g_img_t.*msk;
as(g_img_t)

%% Calculating g-factor Iterative Method
tic

psf = repmat(psf_yz,[1 1 1 nCh]);
% aa = (FFT_1D_Edwin(coilSensNew,'kspace',1));
aa = (FFT_1D_Edwin(CoilSensitivity,'kspace',1));
enc = psf.*aa;
enc = FFT_1D_Edwin(enc,'image',1);
% enc = CoilSensitivity;

% Only collapsed colums and slices...
g_img_i = zeros(N,N,N);
for jj=1:N/Rz  %2D
sli = jj:N/Rz:N;
    for ii=1:N/Ry
        col=ii:N/Ry:N;
        enc_res = single(reshape(permute(enc(:,col,sli,: ),[1 4 2 3]),N*nCh,Ry*Rz));
        for kk=1:N
            enc_res1 = enc(kk,col,sli,:);
            enc_res1 = reshape(permute(enc_res1,[1 4 2 3]),nCh,[]);
            EE = enc_res1'*enc_res1;
%             % Approach 1 Theoretical
%             EE_inv = pinv((EE));
%             EE_diag = diag(EE);
%             EE_inv_diag = abs(diag(EE_inv));
%             g_f = sqrt(EE_diag.*EE_inv_diag);
%             g_f = reshape(g_f,1,Ry,Rz);
%             g_img(kk,col,sli)=g_f;
            % Approach 2 iterative
            g_f = zeros(1,size(EE,1));
            for ll=1:size(EE,1)
                e_rho = zeros(size(EE,1),1);
                e_rho(ll) = 1;
                [d,flag]=lsqr(EE,e_rho);
                g_f(ll) = sqrt(e_rho'*d)*sqrt(e_rho'*EE*e_rho);
            end
                g_f = reshape(g_f,1,Ry,Rz);
                g_img_i(kk,col,sli)=g_f;
        end
    end
end

disp('Iterative method: ')
toc

msk=i_t;
msk(msk~=0) = 1;
g_img_i = g_img_i.*msk;
as(g_img_i)

%% Aproach from internet
% C = squeeze(CoilSensitivity(:,:,40,:));
% RX=1;
% RY=2;
% [NX NY L] = size(C);
% NRX = NX/RX;
% NRY = NY/RY;
% g = zeros(NX,NY);
% nc = zeros(NX,NY);
% for ii=1:NX
%     for jj=1:NY
%         if abs(C(ii,jj,1)) < 1e-6
%             g(ii,jj) = 0;
%         else
%             for LX=0:RX-1
%                 for LY=0:RY-1
%                     ndx = mod((ii-1)+LX*NRX,NX)+1;
%                     ndy = mod((jj-1)+LY*NRY,NY)+1;
%                     CT = C(ndx, ndy, :);
%                     CT = CT(:);
%                     if ((LX==0) && (LY==0))
%                         s = CT;
%                     elseif abs(CT(1)) > 1e-6
%                         s = [s CT];
%                     end
%                 end
%                 nc(ii,jj) = nc(ii,jj)+1;
%             end
%          scs = s'*s;
%         scsi = inv(scs);
%         g(ii,jj) = sqrt(scs(1,1)*scsi(1,1));
%         end
%     end
% end
% as(g)
