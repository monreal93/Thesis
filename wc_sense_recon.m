%% Tasks
%   -  
%   -  How to define the indices of collapes slice group, line 30
%   -  How to define the slice_ind, line 39
%   -  Change code to work with all slices, line 24 only 1 slice
%   - 
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
close all; clear all; clc

load Data/img_wc.mat
load Data/sens_map.mat
load Data/psf_yz.mat
load Data/param.mat

%% SENSE recon for Wave-CAIPI
size_x = size(i_wc, 1)/param.ov;
num_chan = size(i_wc, 4);
y_skip = size(i_wc, 2);

% Selecting slice 1 and CoilSensitv of 3 slices
i_wc = i_wc(:,:,7,:);           
cs = zeros(60,60,3,32);
cs(:,:,1,:) = CoilSensitivity(:,:,7,:);
cs(:,:,2,:) = CoilSensitivity(:,:,20,:);
cs(:,:,3,:) = CoilSensitivity(:,:,60,:);

slice_ind = 7:20:60;      % indices of slices in the collapsed slice group

msk_roi = (cs~=0);
msk_roi = msk_roi(:,:,:,1);

Img_WAVE = zeros(size(i_wc,1)/(param.ov),size(i_wc,2)*(param.Ry),size(i_wc,3)*(param.Rz));

lsqr_iter = 200;
lsqr_tol = 1e-3; 

shift_amount = [-1,0,1] .* size(i_wc,2) / 2;


tic

psf_use = repmat(psf_yz(:,:,slice_ind), [1,1,1,param.num_chan]);
    
receive_use = zeros(size(cs));
mask_use = zeros(size(msk_roi));
    
for nz = 1:param.Rz
    receive_use(:,:,nz,:) = circshift(cs(:,:,nz,:), [0,shift_amount(nz),0,0]);
    psf_use(:,:,nz,:) = circshift(psf_use(:,:,nz,:), [0,shift_amount(nz),0,0]);
    mask_use(:,:,nz) = circshift(msk_roi(:,:,nz), [0,shift_amount(nz),0]);
end
mask_use = sum(mask_use,3);
    
for cey = 1:y_skip

    cey_ind = cey : y_skip : cey + (param.Ry-1) * y_skip;
        
    if sum(sum(mask_use(:,cey_ind),1),2) > 0

        psfs = psf_use(:,cey_ind,:,:);
        rcv = receive_use(:,cey_ind,:,:);

        rhs = squeeze( i_wc(:,cey,:,:) );

        param.psfs = psfs;
        param.rcv = rcv;

        [res, flag, relres, iter] = lsqr(@apply_wave, rhs(:), lsqr_tol, lsqr_iter, [], [], [], param);        

        Img_WAVE(:,cey_ind,:) = reshape(res, [param.img_len, param.Ry, param.Rz]);            
    end
end
    
for nz = 1:param.Rz
    Img_WAVE(:,:,nz) = circshift(Img_WAVE(:,:,nz), [0,-shift_amount(nz),0,0]);
end
    
toc            

as(Img_WAVE)
% mosaic(imrotate(Img_WAVE(30:end-30,7:end-7,:), 90), 1, param.Rz, 3, 'Wave-CAIPI', [0,3e-2])