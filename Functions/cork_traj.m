function [kx,ky,kz]=cork_traj(x,y,z,del_ky,del_kz,i_fov,Ry,Rz,x_points,caipi,caipi_del,plt)

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
kx = zeros(x_points+1,i_fov(2)/Ry,i_fov(3)/Rz);
ky =kx; kz = kx;

y = y-(del_ky/2);
z = z-(del_kz/2);
original_z = z;
kk=1;
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
        if kk<Rz/caipi_del
                z = z - ((del_kz/(Rz/caipi_del))*kk);
%                 z = z - ((del_kz/caipi_del));
                kk=kk+caipi_del;
        else
            kk=1;
        end
      end
end

figure;scatter(kz(:),ky(:)); xlabel('Kz');ylabel('Ky');axis image

%  Ploting all the trajerctories
if plt == 1
 figure;
    for j = 1:i_fov(3)/Rz   
        xx = squeeze(kx(:,:,j));
        yy = squeeze(ky(:,:,j));
        zz = squeeze(kz(:,:,j));
        for l=1:i_fov(2)/Ry
            plot3(xx(:,l),yy(:,l),zz(:,l));
            hold on
        end
    end
    xlabel('Kx');ylabel('Ky');zlabel('Kz');view([90 0 0])
end