function [i_wc_recon,rmse] = wc_sense_recon(i_wc,CoilSensitivity,psf_yz,param)
%% User input
recon_all_sl = true;
% sl_wc = 10;                                % Slice from the undersampled WAVE caipi image, < size of undersampled img.

%% SENSE recon for Wave-CAIPI
y_skip = size(i_wc, 2);
i_wc_recon = complex(zeros(size(CoilSensitivity(:,:,:,1))));
param.psf_len = size(i_wc, 1);
param.img_len = size(i_wc,1)/param.ov;
param.pad_size = (param.psf_len - param.img_len) / 2;

% Reconstructing all phantom or only one slice of i_wc
if recon_all_sl
    last_iter_wc = size(i_wc,3);
else
    last_iter_wc = 1; 
end

for iter_wc = 1:last_iter_wc
    %tic
    if recon_all_sl
        sl_wc = iter_wc;
    end

    i_wc_sl = i_wc(:,:,sl_wc,:); 
    cs = zeros(size(i_wc_sl,1)/param.ov,size(i_wc_sl,2)*param.Ry,param.Rz,param.nCh);
    slice_ind = zeros(1,param.Rz);
    sl_i_t = sl_wc;

    % Extracting CS of only overlaped slices and getting slice indices
    cs(:,:,1,:) = CoilSensitivity(:,:,sl_i_t,:);
    slice_ind(1) = sl_i_t;

    for ii=2:param.Rz
            if (sl_i_t+size(i_wc,3)) > size(i_wc,3)*param.Rz
                sl_i_t = sl_i_t+size(i_wc,3)-size(CoilSensitivity,3);
                cs(:,:,ii,:) = CoilSensitivity(:,:,sl_i_t,:);
                slice_ind(ii)=sl_i_t;
            else
                sl_i_t = sl_i_t+size(i_wc,3);
                cs(:,:,ii,:) = CoilSensitivity(:,:,sl_i_t,:);
                slice_ind(ii)=sl_i_t;
            end
    end
    
        msk_roi = (cs~=0);
        msk_roi = msk_roi(:,:,:,1);
        Img_WAVE = zeros(size(i_wc_sl,1)/(param.ov),size(i_wc_sl,2)*(param.Ry),size(i_wc_sl,3)*(param.Rz));
        lsqr_iter = 200;
        lsqr_tol = 1e-5; 

        shift_amount = zeros(1,param.Rz);
        for ii=1:param.Rz-1
            shift_amount(ii+1) = ii;
        end
        
        if param.N == param.slices
            shift_amount = shift_amount .* (size(i_wc_sl,2)/(param.Ry/param.caipi_del));
        else
            shift_amount = shift_amount .* (size(i_wc_sl,2)/param.Rz);
            shift_amount = floor(shift_amount);
        end
        
        psf_use = repmat(psf_yz(:,:,slice_ind), [1,1,1,param.nCh]);   
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
        
        i_wc_recon(:,:,slice_ind) = Img_WAVE;
        
        %toc
end

%% Calculatin RMSE
i_t=phantom3d(param.N);
msk = CoilSensitivity(:,:,:,1);
msk(msk~=0) = 1;
i_t = double(i_t.*msk);

if param.N ~= param.slices 
    i_t = i_t(:,:,((param.N-param.slices)/2)+1:((param.N-param.slices)/2)+param.slices);
end
i_wc_recon_abs = abs(i_wc_recon);
i_wc_recon_abs = (i_wc_recon_abs-min(i_wc_recon_abs(i_wc_recon_abs~=0)))/(max(i_wc_recon_abs(i_wc_recon_abs~=0))-min(i_wc_recon_abs(i_wc_recon_abs~=0)));

rmse = sqrt(immse(i_wc_recon_abs,i_t)./abs(mean(i_t(:)~=0)));

end
