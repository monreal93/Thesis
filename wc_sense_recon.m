%% Tasks
%   -  
%   -  How to define the indices of collapes slice group, line 30
%   -  How to define the slice_ind, line 39
%   -  Change code to work with all slices, line 24 only 1 slice
%   -  Add logic to wrap slices around if selected sl  does not work...
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
close all; clear all; clc

load Data/img_wc.mat
load Data/sens_map.mat
load Data/psf_yz.mat
load Data/param.mat
load Data/phantom.mat

% i_wc = permute(i_wc,[1 3 2 4]);
% CoilSensitivity = permute(CoilSensitivity,[1 3 2 4]);
CoilSensitivity = CoilSensitivity((((size(i_t,2)*param.ov)-(size(i_t,2)))/2)+1:end-(((size(i_t,2)*param.ov)-(size(i_t,2)))/2),:,:,:);

caipi2d = false;



%% SENSE recon for Wave-CAIPI

size_x = size(i_wc, 1)/param.ov;
num_chan = size(i_wc, 4);
y_skip = size(i_wc, 2);
sl_wc = 22;                                % Slice from the undersampled WAVE caipi image, < size of undersampled img.

% Selecting CoilSensitivity of only overlapped slices
cs = zeros(size(i_wc,1)/param.ov,size(i_wc,2)*param.Ry,param.Rz,param.num_chan);

if param.Rz==1
    sl_i_t = sl_wc;
    cs(:,:,1,:) = CoilSensitivity(:,:,sl_i_t,:);   
else
sl_i_t = (param.Rz*sl_wc)+((size(i_wc,2)/2)-sl_wc);
cs(:,:,1,:) = CoilSensitivity(:,:,sl_i_t,:);
    for ii=2:param.Ry
        if (sl_i_t+size(i_wc,2)) > size(i_wc,2)*param.Ry
            sl_i_t = sl_i_t-size(i_wc,2);
            cs(:,:,ii,:) = CoilSensitivity(:,:,sl_i_t,:);
        else
            sl_i_t = sl_i_t+size(i_wc,2);
            cs(:,:,ii,:) = CoilSensitivity(:,:,sl_i_t,:);
        end
    end
end

i_wc = i_wc(:,:,sl_wc,:); 

slice_ind = sl_wc:size(i_wc,2):size(i_wc,2)*param.Ry;      % indices of slices in the collapsed slice group
% slice_ind = sl*2:size(i_wc,2)/2:(size(i_wc,2)*param.Ry)-1;      % indices of slices in the collapsed slice group
% slice_ind = [7 37];
msk_roi = (cs~=0);
msk_roi = msk_roi(:,:,:,1);

Img_WAVE = zeros(size(i_wc,1)/(param.ov),size(i_wc,2)*(param.Ry),size(i_wc,3)*(param.Rz));

lsqr_iter = 200;
lsqr_tol = 1e-5; 

if param.Rz == 1
    shift_amount = 0;
elseif param.Rz ==2
    shift_amount = [-1,1] .* size(i_wc,2)/ 2;
elseif param.Rz ==3
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

%% SENSE recon for 2D-CAIPI
if caipi2d
    i_wc = i_wc(end/2-29:end/2+30,:,:,:);
    Ry = 2;
    Rz = 2;

    size_x = size(i_wc, 1);
    num_chan = size(i_wc, 4);
    y_skip = size(i_wc, 2);


    i_ind = zeros(num_chan*Ry*Rz*size_x, 1);
    for t = 1:size_x
        ind = 1+(t-1)*num_chan : t*num_chan;
        ind = repmat(ind, [Ry*Rz,1]);
        i_ind(1+(t-1)*Ry*Rz*num_chan : t*Ry*Rz*num_chan) = ind(:);
    end


    j_ind = zeros(num_chan*Ry*Rz*size_x, 1);
    for t = 1:size_x
        ind = 1 + (t-1)*Ry*Rz : t*Ry*Rz;
        ind = repmat(ind(:), [1,num_chan]);    
        j_ind(1+(t-1)*Ry*Rz*num_chan : t*Ry*Rz*num_chan) = ind(:);
    end


    % FOV/2 inter-slice shift
    shift_amount = [-1,0,1] .* size(i_wc, 2) / 2;

    Img_2DCAIPI = zeros(size(msk_roi));
    warning('off','all');


    tic

    receive_us%% SENSE recon for 2D-CAIPI

    i_wc = i_wc(end/2-29:end/2+30,:,:,:);
    Ry = 2;
    Rz = 2;

    size_x = size(i_wc, 1);
    num_chan = size(i_wc, 4);
    y_skip = size(i_wc, 2);


    i_ind = zeros(num_chan*Ry*Rz*size_x, 1);
    for t = 1:size_x
        ind = 1+(t-1)*num_chan : t*num_chan;
        ind = repmat(ind, [Ry*Rz,1]);
        i_ind(1+(t-1)*Ry*Rz*num_chan : t*Ry*Rz*num_chan) = ind(:);
    end


    j_ind = zeros(num_chan*Ry*Rz*size_x, 1);
    for t = 1:size_x
        ind = 1 + (t-1)*Ry*Rz : t*Ry*Rz;
        ind = repmat(ind(:), [1,num_chan]);    
        j_ind(1+(t-1)*Ry*Rz*num_chan : t*Ry*Rz*num_chan) = ind(:);
    end


    % FOV/2 inter-slice shift
    shift_amount = [-1,0,1] .* size(i_wc, 2) / 2;

    Img_2DCAIPI = zeros(size(msk_roi));
    warning('off','all');


    tic

    receive_use = zeros(size(cs));
    mask_use = zeros(size(msk_roi));

    for nz = 1:Rz
        receive_use(:,:,nz,:) = circshift(cs(:,:,nz,:), [0,shift_amount(nz),0,0]);
        mask_use(:,:,nz) = circshift(msk_roi(:,:,nz,:), [0,shift_amount(nz),0,0]);
    end
    mask_use = sum(mask_use,3);

    for cey = 1:y_skip
        % collapsing indices
        cey_ind = cey : y_skip : cey + (Ry-1) * y_skip;

        if sum(sum(mask_use(:,cey_ind),1),2) > 0

            data = permute(squeeze(i_wc(:,cey,:,:)), [2,1]);
            Data = double(data(:));     % (coil*x_size) X 1

            encoding = receive_use(:,cey_ind,:,:);
            encoding = double(permute(encoding, [2,3,4,1]));

            Encoding = sparse(i_ind, j_ind, encoding(:), num_chan*size_x, Ry*Rz*size_x, num_chan*Ry*Rz*size_x);

            res = Encoding \ Data;
            res = permute(reshape(res, [Ry, Rz, size_x]), [3,1,2]);
            Img_2DCAIPI(:,cey_ind,:) = res;

        end
    end

    for nz = 1:Rz
        % undo slice shift
        Img_2DCAIPI(:,:,nz) = circshift(Img_2DCAIPI(:,:,nz), [0,-shift_amount(nz),0,0]);
    end
    toc            

    as(Img_2DCAIPI) = zeros(size(cs));
    mask_use = zeros(size(msk_roi));

    for nz = 1:Rz
        receive_use(:,:,nz,:) = circshift(cs(:,:,nz,:), [0,shift_amount(nz),0,0]);
        mask_use(:,:,nz) = circshift(msk_roi(:,:,nz,:), [0,shift_amount(nz),0,0]);
    end
    mask_use = sum(mask_use,3);

    for cey = 1:y_skip
        % collapsing indices
        cey_ind = cey : y_skip : cey + (Ry-1) * y_skip;

        if sum(sum(mask_use(:,cey_ind),1),2) > 0

            data = permute(squeeze(i_wc(:,cey,:,:)), [2,1]);
            Data = double(data(:));     % (coil*x_size) X 1

            encoding = receive_use(:,cey_ind,:,:);
            encoding = double(permute(encoding, [2,3,4,1]));

            Encoding = sparse(i_ind, j_ind, encoding(:), num_chan*size_x, Ry*Rz*size_x, num_chan*Ry*Rz*size_x);

            res = Encoding \ Data;
            res = permute(reshape(res, [Ry, Rz, size_x]), [3,1,2]);
            Img_2DCAIPI(:,cey_ind,:) = res;

        end
    end

    for nz = 1:Rz
        % undo slice shift
        Img_2DCAIPI(:,:,nz) = circshift(Img_2DCAIPI(:,:,nz), [0,-shift_amount(nz),0,0]);
    end
    toc            

    as(Img_2DCAIPI)
end