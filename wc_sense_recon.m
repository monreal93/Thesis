%% Tasks
%   -  
%   -  
%   -  
%   - 
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
%close all; clear all; clc
% load /home/amonreal/Documents/Thesis/GitHub/Thesis/Data/img_wc.mat
% load /home/amonreal/Documents/Thesis/GitHub/Thesis/Data/sens_map.mat
% load /home/amonreal/Documents/Thesis/GitHub/Thesis/Data/psf_yz.mat
% load /home/amonreal/Documents/Thesis/GitHub/Thesis/Data/param.mat
% load /home/amonreal/Documents/Thesis/GitHub/Thesis/Data/phantom.mat

%% User input
caipi2d = false;
recon_all_sl = false;
sl_wc = 7;                                % Slice from the undersampled WAVE caipi image, < size of undersampled img.

%% SENSE recon for Wave-CAIPI
size_x = size(i_wc, 1)/param.ov;
num_chan = size(i_wc, 4);
y_skip = size(i_wc, 2);
OutputImg = complex(zeros(size(CoilSensitivity(:,:,:,1))));

% Reconstructing all phantom or only one slice of i_wc
if recon_all_sl
    last_iter_wc = size(i_wc,3);
else
    last_iter_wc = 1; 
end

for iter_wc = 1:last_iter_wc
    if recon_all_sl
        sl_wc = iter_wc;
    end

    i_wc_sl = i_wc(:,:,sl_wc,:); 
    cs = zeros(size(i_wc_sl,1)/param.ov,size(i_wc_sl,2)*param.Ry,param.Rz,param.num_chan);
    slice_ind = zeros(1,param.Rz);

    % Checking starting slice on full img, depending on R odd/pair
    if  mod(param.Ry,2)
        % R is odd
        sl_i_t=sl_wc;
    else
        % R is pair
        sl_i_t = sl_wc+(size(i_wc,2)/2);
    end

    % Extracting CS of only overlaped slices and getting slice indices
    cs(:,:,1,:) = CoilSensitivity(:,:,sl_i_t,:);
    slice_ind(1) = sl_i_t;
    for ii=2:param.Ry
            if (sl_i_t+size(i_wc_sl,2)) > size(i_wc_sl,2)*param.Ry
                sl_i_t = sl_i_t-size(i_wc_sl,2);
                cs(:,:,ii,:) = CoilSensitivity(:,:,sl_i_t,:);
                slice_ind(ii)=sl_i_t;
            else
                sl_i_t = sl_i_t+size(i_wc_sl,2);
                cs(:,:,ii,:) = CoilSensitivity(:,:,sl_i_t,:);
                slice_ind(ii)=sl_i_t;
            end
    end

        msk_roi = (cs~=0);
        msk_roi = msk_roi(:,:,:,1);
        Img_WAVE = zeros(size(i_wc_sl,1)/(param.ov),size(i_wc_sl,2)*(param.Ry),size(i_wc_sl,3)*(param.Rz));
        lsqr_iter = 200;
        lsqr_tol = 1e-5; 

        if param.Rz == 1
            shift_amount = 0;
        elseif param.Rz ==2
            shift_amount = [-1,0] .* size(i_wc_sl,2)/ 2;
        elseif param.Rz ==3
            shift_amount = [-1,0,1] .* size(i_wc_sl,2)/ 2;
        elseif param.Rz == 4
            shift_amount = [-2,-1,0,1] .* (size(i_wc_sl,2)/ 2);
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

        for cey = 1:y_skip
            cey_ind = cey : y_skip : cey + (param.Ry-1) * y_skip;
            if sum(sum(mask_use(:,cey_ind),1),2) > 0
                psfs = psf_use(:,cey_ind,:,:);
                rcv = receive_use(:,cey_ind,:,:);
                rhs = squeeze( i_wc_sl(:,cey,:,:) );
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
        OutputImg(:,:,slice_ind) = Img_WAVE;    
end

% Plotting all volume or only some slices
if recon_all_sl
    as(OutputImg)
else
    as(Img_WAVE)
end

%% Calculatin RMSE
aa=phantom3d(N);
bb=aa(:,:,slice_ind);
rmse=immse(Img_WAVE,bb);
fprintf('\n The mean-squared error is %0.4f\n', rmse);

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

    as(Img_2DCAIPI)
end



