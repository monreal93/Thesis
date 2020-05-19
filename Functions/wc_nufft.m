function kspace_nufft = wc_nufft(i_t,kx,ky,kz,param)

% Normalize trajectories
tr = zeros(3,length(kx(:)));
tr(1,:) = kx(:)./param.k_fov(1);
tr(2,:) = ky(:)./param.k_fov(2);
tr(3,:) = kz(:)./param.k_fov(3);

osf = 2; wg = 3; sw = 8;

FT = gpuNUFFT(tr,col(ones(size(kx(:)))),osf,wg,sw,[param.N,param.N,param.slices],[],true);
kspace_nufft = FT*i_t;
if param.plt == 1
    figure; view(2)
    scatter3(kx(:),ky(:),kz(:),ones(size(kz(:))).*50,abs(kspace_nufft(:,1)),'.')
    xlabel('Kx');ylabel('Ky');zlabel('Kz'); view([90 0 0])
end
clear ky kz

kspace_nufft = reshape(kspace_nufft,param.N*param.ov,param.N/param.Ry,param.slices/param.Rz,param.nCh);
kspace_nufft = flip(kspace_nufft,2);
kspace_nufft = flip(kspace_nufft,3);

% Adding noise to k-space data
kspace_nufft = addGaussianNoise_cmplx(kspace_nufft,20);





