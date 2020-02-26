%% Tasks
%   -  
%   -  How to define the indices of collapes slice group, line 30
%   -  How to define the slice_ind, line 39
%   -  Change code to work with all slices, line 24 only 1 slice
%   -  Add logic to wrap slices around if selected sl  does not work...
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
% close all; clear all; clc

load Data/img_wc.mat
% load Data/sens_map.mat
% load Data/psf_yz.mat
% load Data/param.mat
% load Data/phantom.mat
% aa = i_wc; i_wc = aa;

%% SENSE recon for Wave-CAIPI
size_x = size(i_wc, 1)/param.ov;
num_chan = size(i_wc, 4);
y_skip = size(i_wc, 2);
sl_wc = 22;                                % Slice from the undersampled WAVE caipi image, < size of undersampled img.

% Selecting CoilSensitivity of only overlapped slices
cs = zeros(size(i_wc,2)*param.Ry,size(i_wc,3)*param.Rz,Ry,param.num_chan);
sl_i_t = (param.Rz*sl_wc)+((size(i_wc,2)/2)-sl_wc);
cs(:,:,1,:) = CoilSensitivity(:,:,sl_i_t,:);
if param.Ry>1
    for ii=2:param.Ry
        if (sl_i_t+size(i_wc,2)) > size(i_wc,2)*param.Ry
            cs(:,:,ii,:) = CoilSensitivity(:,:,sl_i_t-size(i_wc,2),:);
        else
            cs(:,:,ii,:) = CoilSensitivity(:,:,sl_i_t+size(i_wc,2),:);
        end
    end
end

i_wc = i_wc(:,:,sl_wc,:); 

slice_ind = sl_wc:size(i_wc,2):size(i_wc,2)*param.Ry;      % indices of slices in the collapsed slice group
% slice_ind = sl*2:size(i_wc,2)/2:(size(i_wc,2)*param.Ry)-1;      % indices of slices in the collapsed slice group
% slice_ind = [30 45];
msk_roi = (cs~=0);
msk_roi = msk_roi(:,:,:,1);

Img_WAVE = zeros(size(i_wc,1)/(param.ov),size(i_wc,2)*(param.Ry),size(i_wc,3)*(param.Rz));

lsqr_iter = 200;
lsqr_tol = 1e-3; 

if param.Ry == 1
    shift_amount = 0;
elseif param.Ry ==2
    shift_amount = [-1,1] .* size(i_wc,2)/ 2;
elseif param.Ry ==3
    shift_amount = [-1,0,1] .* size(i_wc,2)/ 2;
end

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
% mask_use = circshift(mask_use,-30,2);
    
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