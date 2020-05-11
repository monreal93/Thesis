function g_av_t=gfactorcalc(x,p_bw,p_s,k_fov,i_fov,Ry,Rz,x_points,gamma,z_offset,coil_radius,loop_radius,Resolution,i_t,nCh,plt,N,slices,ov,gf)
  
%% Calculating parameters for Gradients
[~,~,~,x_calc,y_calc,z_calc,~,~,t,t_prime]=grad_param(x(1),x(2),x(3),x(4),p_bw,p_s,k_fov,i_fov,Ry,Rz,x_points,gamma);

%% Generate coil sensitivity
CoilSensitivity = coil_sens(z_offset,coil_radius,loop_radius,Resolution,i_t,nCh,plt);

%% Cropping image and CoilSens, if z <> x,y
i_t = i_t(:,:,((N-slices)/2)+1:((N-slices)/2)+slices);
CoilSensitivity = CoilSensitivity(:,:,((N-slices)/2)+1:((N-slices)/2)+slices,:);

%% Creating PSFs
[~,~,psf_yz]=psf(N,slices,i_fov,k_fov,ov,p_s,x_points,i_t,t,t_prime,x_calc,y_calc,z_calc,plt);

%% Calculating G-Factor
if gf == 1 || gf==3
    [~,g_av_t,g_max_t]=gfact_teo(N,slices,Ry,Rz,nCh,i_t,psf_yz,CoilSensitivity);
end
if gf == 2 || gf==3
    [~,g_av_i,g_max_i]=gfact_iter(N,Ry,Rz,nCh,i_t,psf_yz,CoilSensitivity);
end
end