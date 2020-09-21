% Script to generate the Wave-CAIPI image using the NUFFT function

function kspace_nufft = wc_nufft(i_t,kx,ky,kz,param)

% Normalize trajectories
tr = zeros(3,length(kx(:)));
tr(1,:) = kx(:)./param.k_fov(1);
tr(2,:) = ky(:)./param.k_fov(2);
tr(3,:) = kz(:)./param.k_fov(3);

if param.gpu
    osf = 2; wg = 3; sw = 8;
    FT = gpuNUFFT(tr,col(ones(size(kx(:)))),osf,wg,sw,[param.Nx,param.Ny,param.slices],[],true);
    kspace_nufft1=FT*i_t;
else
    w=1;
    FT = NUFFT3D(tr', col(sqrt(w)), 1, 0, [param.Nx,param.Ny,param.slices], 2);
    kspace_nufft1 = zeros(param.Nx*param.Ny*param.slices*param.ov/(param.Ry*param.Rz),param.nCh);
    for ii=1:param.nCh
        kspace_nufft1(:,ii) = FT*i_t(:,:,:,ii);
    end
end

if param.plt == 1
    figure; view(2)
    scatter3(kx(:),ky(:),kz(:),ones(size(kz(:))).*50,abs(kspace_nufft1(:,1)),'.')
    xlabel('Kx');ylabel('Ky');zlabel('Kz'); view([90 0 0])
end
clear ky kz

kspace_nufft = reshape(kspace_nufft1,param.Nx*param.ov,param.Ny/param.Ry,param.slices/param.Rz,param.nCh);
kspace_nufft = flip(kspace_nufft,2);
kspace_nufft = flip(kspace_nufft,3);

% Adding noise to k-space data
kspace_nufft = addGaussianNoise_cmplx(kspace_nufft,20);





