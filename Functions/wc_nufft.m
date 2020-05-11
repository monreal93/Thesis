function kspace_nufft = wc_nufft(kx,ky,kz,k_fov,N,slices,i_t,ov,Ry,Rz,nCh,plt)

% Normalize trajectories
tr = zeros(3,length(kx(:)));
tr(1,:) = kx(:)./k_fov(1);
tr(2,:) = ky(:)./k_fov(2);
tr(3,:) = kz(:)./k_fov(3);

%clear ky kz
% %%% Density compensation
% tr_notScaled = zeros(3,length(kx(:)));
% tr_notScaled (1,:) = kx(:)./1;
% tr_notScaled (2,:) = ky(:)./1;
% tr_notScaled (3,:) = kz(:)./1;
% w = iterative_dcf_estimation(tr_notScaled,6,2.1); 

disp('Generate NUFFT Operator without coil sensitivities');
disp_slice=N/2;
useGPU = true;
useMultiCoil = 0;
osf = 2; wg = 3; sw = 8;
imwidth = N;

% % Apprach 1:
% i_t = padarray(i_t,(((size(i_t,2)*ov)-size(i_t,2))/2),'both');
% FT = gpuNUFFT(tr,col(ones(size(kx(:)))),osf,wg,sw,[N*ov,N,N],[],true);
% % FT = gpuNUFFT(tr,w,osf,wg,sw,[N*ov,N,N],[],true);
% kspace_nufft = FT*i_t;
% if plt == 1
%     figure; view(2)
%     scatter3(kx(:),ky(:),kz(:),ones(size(kz(:))).*50,abs(kspace_nufft(:,1)),'.')
%     xlabel('Kx');ylabel('Ky');zlabel('Kz'); view([90 0 0])
% end
% kspace_nufft = reshape(kspace_nufft,N*ov,N/Ry,N/Rz,nCh);
% kspace_nufft = flip(kspace_nufft,2);
% kspace_nufft = flip(kspace_nufft,3);
% % i_wc = imag_wc(N,Ry,Rz,caipi,kspace_nufft,plt);

% Approach 2:
FT = gpuNUFFT(tr,col(ones(size(kx(:)))),osf,wg,sw,[N,N,slices],[],true);
kspace_nufft = FT*i_t;
if plt == 1
    figure; view(2)
    scatter3(kx(:),ky(:),kz(:),ones(size(kz(:))).*50,abs(kspace_nufft(:,1)),'.')
    xlabel('Kx');ylabel('Ky');zlabel('Kz'); view([90 0 0])
end
clear ky kz

kspace_nufft = reshape(kspace_nufft,N*ov,N/Ry,slices/Rz,nCh);
kspace_nufft = flip(kspace_nufft,2);
kspace_nufft = flip(kspace_nufft,3);
% i_wc1 = imag_wc(N,Ry,Rz,caipi,kspace_nufft1,plt);







