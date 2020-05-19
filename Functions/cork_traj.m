function [kx,ky,kz,r_ps]=cork_traj(x,y,z,del_ky,del_kz,param)

%% Moving in y direction first and then z...
% kx = zeros(x_points+1,i_fov(2)/Ry,i_fov(3)/Rz);
% ky =kx; kz = kx;
% 
% y = y-(del_ky/2);
% z = z-(del_kz/2);
% % original_z = z;
% original_y = y;
% kk=1;
% for j=1:(i_fov(2)/Ry)
%     for l=1:(i_fov(3)/Rz)
%             kx(:,j,l) = x;
%             ky(:,j,l) = y;
%             kz(:,j,l) = z;
% %             z = z -(del_kz);
%             y=y-(del_ky);
%     end
% %           y = y-(del_ky);
% %           z = original_z;
%           z = z-(del_kz);
%           y = original_y;
%       if caipi
%         if kk<Ry/caipi_del
%                 y = y - ((del_ky/(Ry/caipi_del))*kk);
%                 kk=kk+caipi_del;
%         else
%             kk=1;
%         end
%       end
% end

%% Different approach moving in z direction first and then y....
% Creating null matrices for K points
kx = zeros(param.x_points+1,param.i_fov(2)/param.Ry,param.i_fov(3)/param.Rz);
ky =kx; kz = kx;

y = y-(del_ky/2);
z = z-(del_kz/2);
original_z = z;
kk=1;
for j=1:(param.i_fov(2)/param.Ry)
    for l=1:(param.i_fov(3)/param.Rz)
            kx(:,j,l) = x;
            ky(:,j,l) = y;
            kz(:,j,l) = z;
            z = z -(del_kz);
    end
          y = y-(del_ky);
          z = original_z;
      if param.caipi
        if kk<param.Rz/param.caipi_del
                z = z - ((del_kz/(param.Rz/param.caipi_del))*kk);
                kk=kk+param.caipi_del;
        else
            kk=1;
        end
      end
end

ky=ky-mean(ky(:));
kz=kz-mean(kz(:));

% figure;scatter(kz(:),ky(:),'.');title('k-space trajectory');xlabel('Kz');ylabel('Ky');axis image

% Calculating real pixel size
r_ps = [1./(max(kx(:))-min(kx(:))) 1./(max(ky(:))-min(ky(:))) 1./(max(kz(:))-min(kz(:)))];

%  Ploting all the trajerctories
if param.plt == 1
 figure;
    for j = 1:param.i_fov(3)/param.Rz   
        xx = squeeze(kx(:,:,j));
        yy = squeeze(ky(:,:,j));
        zz = squeeze(kz(:,:,j));
        for l=1:param.i_fov(2)/param.Ry
            plot3(xx(:,l),yy(:,l),zz(:,l));
            hold on
        end
    end
    xlabel('Kx');ylabel('Ky');zlabel('Kz');view([90 0 0])
end