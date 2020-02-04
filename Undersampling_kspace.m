addpath(genpath('/home/amonreal/Documents/Thesis/Matlab_scripts/Code4Alejandro/shepp_logan3d'))

%%% Code to use when a NIFTI is used.
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
%%%%%%%%%%%%%%%%%%%%%%%

%%% Using a Phantom....
i_t = phantom3d(240);
i_f = fftshift(fftn(i_t));
i_fov = [size(i_t,1) size(i_t,2) size(i_t,3)];                                   % mm3
p_s = [1e-3 1e-3 1e-3];
voxel_size = p_s(1)*p_s(2)*p_s(3)*.001;
x_points = i_fov(1)*6;
gy_gz_amplit = 6e-3;
sins = 7;
p_bw = 70;
Ry = 3; Rz = 3;
k_fov = 1./p_s;
%%%%%%%%%%%%%

t_r = 1/p_bw;
delta_t = 1e-7;
t_prime = delta_t:delta_t:t_r;
t = 0+(t_r/x_points):t_r/x_points:t_r;
f = sins/t_r;
caipi = 1;                % 1 to dephase pair lines as in 2D CAIPI
gamma = 42.58e6;

del_ky = k_fov(2)/(i_fov(2)/Ry);
del_kz = k_fov(3)/(i_fov(3)/Rz);
t_r = 1/p_bw;
omg = 2*pi*f;

kx = zeros(round(k_fov(2)/del_ky),round(k_fov(3)/del_kz),x_points);
ky =kx; kz = kx;

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

del_ky = (k_fov(2)/(i_fov(2)/Ry))+(max(y)-min(y));
del_kz = (k_fov(3)/(i_fov(3)/Rz))+(max(z)-min(z));
original_y = y; 

for j=1:(k_fov(3)/del_kz)
    for i=1:(k_fov(2)/del_ky)
            kx(i,j,:) = x;
            ky(i,j,:) = y;
            kz(i,j,:) = z;
            y = y -(del_ky);
%         if (caipi == 1) &&(rem(j,2) ==0)&& (((k_fov/del_ky)-i)==1)
%            y = y -del_ky;
%             break 
%         end
    end
          z = z-(del_kz);
          y = original_y;
      if caipi == 1
         if (rem(j,2) ==1)
                y = y - (del_ky/2);
         end
      end
end

%%% ploting all the trajerctories
for j = 1:k_fov/del_ky    
    xx = squeeze(kx(j,:,:));
    yy = squeeze(ky(j,:,:));
    zz = squeeze(kz(j,:,:));
    for i=1:k_fov/del_kz
        plot3(xx(i,:),yy(i,:),zz(i,:));
        hold on
    end
end
xlabel('Kx');ylabel('Ky');zlabel('Kz');view(3)

% %%% Interpolating
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

