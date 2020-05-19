function [psf_y,psf_z,psf_yz]=psf(i_t,t,t_prime,x_calc,y_calc,z_calc,param)

%% Creating cartersian trajectory
del_c = 1./(param.i_fov.*param.p_s);
x_c = zeros(param.x_points+1,size(i_t,2),size(i_t,3));
y_c = x_c;  z_c = x_c;

xc = (param.k_fov(1)/2*-1)+((del_c(1)/param.ov)/2):del_c(1)/param.ov:param.k_fov(1)/2;
yc = repmat((param.i_fov(2)/2*-1),1,size(xc,2))*param.p_s(2)+(param.p_s(2)/2);
zc = repmat((param.i_fov(3)/2*-1),1,size(xc,2))*param.p_s(3)+(param.p_s(3)/2);
original_yc = yc;

for j=1:param.i_fov(3)
    for l=1:param.i_fov(2)
            x_c(:,l,j) = xc;
            y_c(:,l,j) = yc;
            z_c(:,l,j) = zc;
            yc = yc +param.p_s(2);
    end
          zc = zc+param.p_s(3);
          yc = original_yc;
end

% Ploting all the trajerctories
if param.plt == 1
figure;
    for j = 1:param.i_fov(3)   
        for l=1:param.i_fov(2)
            xx = x_c(:,l,j); xx = xx(:)';
            yy = y_c(:,l,j); yy = yy(:)';
            zz = z_c(:,l,j); zz = zz(:)';
            plot3(xx,yy,zz);
            hold on
        end
    end
    xlabel('Kx');ylabel('y');zlabel('z');view(3); title('Hybrid K-space');
end

%% Creating PSF
psf_y = zeros(param.x_points+1,param.N);
psf_z = zeros(param.x_points+1,param.slices);

Px = interp1(t_prime,x_calc,t,'linear','extrap');
Py = interp1(t_prime,y_calc,t,'linear','extrap');
Pz =  interp1(t_prime,z_calc,t,'linear','extrap');

if param.plt==1
    figure; plot(Px,Py)
    hold on
    plot(Px,Pz)
end

% psf Y
    for l=1:param.N
        yy = y_c(:,l,1); yy = yy(:)';
        psf_y(:,l) = exp(-1i*2*pi*Py.*yy);
    end
    
% psf Z
    for l=1:param.slices
        zz = z_c(:,1,l); zz = zz(:)';
        psf_z(:,l) = exp(-1i*2*pi*Pz.*zz);
    end 

 psf_yz = repmat(psf_y,[1,1,size(psf_z,2)]) .* repmat(permute(psf_z, [1,3,2]), [1,size(psf_y,2),1]);
 

%%% Analyzing the power spectrum
psf_f=FFT_1D_Edwin(psf_yz,'image',1);
psf_f1 = psf_f(:,:,1);
t_r = 1/param.p_bw;
t_smp = param.ov*param.N;
fq = t_smp./t_r;
x=psf_f1(:,1);
x=FFT_2D_Edwin(x,'kspace',2,3);
x=x(:); x = x-mean(x);
n=param.N*param.ov;
f=(-n/2:n/2-1)*(fq/n);
power=abs(x).^2/n;
figure;subplot(2,1,1); imagesc(abs(psf_f1.')); title('Psf')
sz=param.N*param.p_bw;
subplot(2,1,2); plot(f,power);  title('Power spectrum'); xlabel('Frequency Hz');

%%%

if param.plt == 1
    figure; imagesc(angle(psf_y).');colormap jet, axis image off ; title('PSF Y')
    figure; imagesc(angle(psf_z).');colormap jet, axis image off; title('PSF Z')
end