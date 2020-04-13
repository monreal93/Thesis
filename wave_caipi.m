

load img_2dcaipi_Ry3_Rz3
%% SENSE recon for 2D-CAIPI

Ry = 3;
Rz = 3;

size_x = size(img_2dcaipi, 1);
num_chan = size(img_2dcaipi, 4);
y_skip = size(img_2dcaipi, 2);

%% load coil sensitivities and roi mask

load receive
load mask_roi

%% load Point Spread Function of Gy and Gz gradients (PSF-Y and PSF-Z) 

os_size = 580;

load PsfY_fit
load PsfZ_fit

figure(11), imagesc( angle(PsfY_fit).' ), colormap jet, axis image off; title('PSF Y')
figure(12), imagesc( angle(PsfZ_fit).' ), colormap jet, axis image off; title('PSF Z')

%% Wave-CAIPI : Ry=3 & Rz=3
%% load undersampled slice group acquired with Wave-CAIPI

PsfYZ_fit = repmat(PsfY_fit,[1,1,size(PsfZ_fit,2)]) .* repmat(permute(PsfZ_fit, [1,3,2]), [1,size(PsfY_fit,2),1]);

load img_wave_Ry3_Rz3


%% SENSE recon for Wave-CAIPI

param = [];
param.psf_len = size(img_wave, 1);                              % psf length
param.img_len = size(mask_roi, 1);                             % image length (without oversampling)
param.pad_size = (param.psf_len - param.img_len) / 2;          % zero pad size (psf length - image length) / 2
param.num_chan = num_chan;
param.Ry = Ry;
param.Rz = Rz;

Img_WAVE = zeros(size(mask_roi));

lsqr_iter = 200;
lsqr_tol = 1e-3; 

shift_amount = [-1,0,1] .* size(img_wave,2) / 2;
slice_ind = 14:40:120;      % indices of slices in the collapsed slice group

tic

psf_use = repmat(PsfYZ_fit(:,:,slice_ind), [1,1,1,num_chan]);
    
receive_use = zeros(size(receive));
mask_use = zeros(size(mask_roi));
    
for nz = 1:Rz
    receive_use(:,:,nz,:) = circshift(receive(:,:,nz,:), [0,shift_amount(nz),0,0]);
    psf_use(:,:,nz,:) = circshift(psf_use(:,:,nz,:), [0,shift_amount(nz),0,0]);
    mask_use(:,:,nz) = circshift(mask_roi(:,:,nz), [0,shift_amount(nz),0]);
end
mask_use = sum(mask_use,3);
    
for cey = 1:y_skip

    cey_ind = cey : y_skip : cey + (Ry-1) * y_skip;
        
    if sum(sum(mask_use(:,cey_ind),1),2) > 0

        psfs = psf_use(:,cey_ind,:,:);
        rcv = receive_use(:,cey_ind,:,:);

        rhs = squeeze( img_wave(:,cey,:,:) );

        param.psfs = psfs;
        param.rcv = rcv;

        [res, flag, relres, iter] = lsqr(@apply_wave, rhs(:), lsqr_tol, lsqr_iter, [], [], [], param);        

        Img_WAVE(:,cey_ind,:) = reshape(res, [param.img_len, Ry, Rz]);            
    end
end
    
for nz = 1:Rz
    Img_WAVE(:,:,nz) = circshift(Img_WAVE(:,:,nz), [0,-shift_amount(nz),0,0]);
end
    
toc            

mosaic(imrotate(Img_WAVE(30:end-30,7:end-7,:), 90), 1, Rz, 3, 'Wave-CAIPI', [0,3e-2])