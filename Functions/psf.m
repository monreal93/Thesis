function [psf_y,psf_z,psf_yz]=psf(N,slices,i_fov,k_fov,ov,p_s,x_points,i_t,t,t_prime,x_calc,y_calc,z_calc,plt)

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
    for j = 1:i_fov(3)   
        for l=1:i_fov(2)
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
% psf_y = zeros(x_points+1,size(ky,2)*Ry);
% psf_z = zeros(x_points+1,size(kz,3)*Rz);

psf_y = zeros(x_points+1,N);
psf_z = zeros(x_points+1,slices);

Px = interp1(t_prime,x_calc,t,'linear','extrap');
Py = interp1(t_prime,y_calc,t,'linear','extrap');
Pz =  interp1(t_prime,z_calc,t,'linear','extrap');

if plt==1
    figure; plot(Px,Py)
    hold on
    plot(Px,Pz)
end

% psf Y
    for l=1:N
        yy = y_c(:,l,1); yy = yy(:)';
        psf_y(:,l) = exp(-1i*2*pi*Py.*yy);
    end
    
% psf Z
    for l=1:slices
        zz = z_c(:,1,l); zz = zz(:)';
        psf_z(:,l) = exp(-1i*2*pi*Pz.*zz);
    end 

 psf_yz = repmat(psf_y,[1,1,size(psf_z,2)]) .* repmat(permute(psf_z, [1,3,2]), [1,size(psf_y,2),1]);
 
 
if plt == 1
    figure; imagesc(angle(psf_y).');colormap jet, axis image off ; title('PSF Y')
    figure; imagesc(angle(psf_z).');colormap jet, axis image off; title('PSF Z')
end