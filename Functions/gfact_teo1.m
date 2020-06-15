function [g_img_t,g_av_t,g_max_t]=gfact_teo1(i_t,psf_yz,CoilSensitivity,param)

        psf = repmat(psf_yz,[1 1 1 param.nCh]);
        CoilSensitivity = padarray(CoilSensitivity,((param.N*param.ov)-param.N)/2,'both');
        
        %%% Testing:
%         psf = conj(psf);
%         CoilSensitivity = conj(CoilSensitivity);
        %%%
        aa = (FFT_1D_Edwin(CoilSensitivity,'kspace',1));
        enc = psf.*aa;
        enc = FFT_1D_Edwin(enc,'image',1);
        
        enc = enc((((param.N*param.ov)-param.N)/2)+1:(((param.N*param.ov)-param.N)/2)+param.N,:,:,:);
        
% Only collapsed colums and slices, pixel by pixel
g_img_t = zeros(param.N,param.N,param.slices);
col = zeros(floor(param.Rz/param.caipi_del),param.Ry);
enc_res = zeros(param.N*param.nCh,param.Ry*param.Rz);

for mm=1:param.Rz-1
    shift_amount(mm+1) = mm;
end
shift_amount = shift_amount .* ((param.N/param.Rz)/(param.Rz/param.caipi_del));
shift_amount = floor(shift_amount);

   for jj=1:param.slices/param.Rz 
        sli = jj:param.slices/param.Rz:param.slices;
        sli = reshape(sli,param.Rz,[]);
        for ii=1:param.N/param.Ry
                        
                for mm=1:param.Rz/param.caipi_del
                    aa = ii+shift_amount(mm):param.N/param.Ry:param.N+(param.N/param.Rz);
                    aa = aa(1:param.Ry);
                    col(mm,:) = aa;
                    if max(col(:)) > param.N
                        idx = find(col>param.N);
                        col(idx) = col(idx)-param.N;
                    end
                end
                col = repmat(col,param.caipi_del,1);
                for mm=1:param.Rz
                    enc_res(:,mm:param.Rz:end) = double(reshape(permute(enc(:,col(mm,:),sli(mm,:),:),[1 4 2 3]),param.N*param.nCh,param.Ry));
                end
                for kk=1:param.N
                        enc_res1 = enc_res(kk:param.N:end,:);
                        EE = enc_res1'*enc_res1;
                        % Approach 1 Theoretical
                        EE_inv = pinv((EE));
                        EE_diag = diag(EE);
                        EE_inv_diag = (diag(EE_inv));
                        g_f = sqrt(EE_diag.*EE_inv_diag);
                        g_f = reshape(g_f,param.Rz,param.Ry);
                            for mm=1:param.Rz
                                g_img_t(kk,col(mm,:),sli(mm,:))=g_f(mm,:);
                            end
                end
                col = zeros(floor(param.Rz/param.caipi_del),param.Ry);
        end
    end
       

        msk=i_t(:,:,:,1);
        msk(msk~=0) = 1;
        g_img_t = g_img_t.*msk;
        g_av_t = abs(mean(g_img_t(g_img_t~=0)));
        g_max_t = abs(max(g_img_t,[],'all'));