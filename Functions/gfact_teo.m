function [g_img_t,g_av_t,g_max_t]=gfact_teo(i_t,psf_yz,CoilSensitivity,param)

        %tic
        
        psf = repmat(psf_yz,[1 1 1 param.nCh]);
        CoilSensitivity = padarray(CoilSensitivity,((param.N*param.ov)-param.N)/2,'both');
        aa = (FFT_1D_Edwin(CoilSensitivity,'kspace',1));
        enc = psf.*aa;
        enc = FFT_1D_Edwin(enc,'image',1);
        % 2D CAIPI effect
%         enc_test1 = enc;
%         enc_test1 = circshift(enc_test1,N/Rz/2,3);
%         for ii=1:param.slices/param.Rz
%             if param.Rz==2
%                 shift=param.slices-(param.slices/(param.Rz*param.Rz));
%             elseif param.Rz==4
%                 shift=(param.slices/param.Rz)*2;
%             elseif param.Rz==1
%                 shift=1;
%             elseif param.Rz==3
%                 F
%             end
%             enc_test1(:,:,ii:shift:param.slices,:) = circshift(enc_test1(:,:,ii:shift:param.slices,:),param.N/2,2);
%         end
%         enc = enc_test1;
        
% Only collapsed colums and slices, pixel by pixel
g_img_t = zeros(param.N,param.N,param.slices);
col = zeros(param.Ry/param.caipi_del,param.Rz);
enc_res = zeros(param.N*param.ov*param.nCh,param.Ry*param.Rz);

for mm=1:param.Rz-1
    shift_amount(mm+1) = mm;
end
shift_amount = shift_amount .* ((param.N/param.Rz)/(param.Rz/param.caipi_del));
            
   for jj=1:param.slices/param.Rz 
        sli = jj:param.slices/param.Rz:param.slices;
        sli = reshape(sli,param.Rz,[]);
        for ii=1:param.N/param.Ry
                        
%                         col1=ii+(param.N/param.Ry/2):param.N/param.Ry:param.N+(param.N/param.Ry/2);
%                         col2=ii:param.N/param.Ry:param.N; 
%                         if max(col1) > param.N
%                             col1(col1>param.N) = (col1(col1>param.N))-param.N;
%                         end
%                     
%                         enc_res_1 = double(reshape(permute(enc(:,col2,sli(1:param.Rz/2:end-1),: ),[1 4 2 3]),param.N*param.ov*param.nCh,param.Rz*param.Ry/2));
%                         enc_res_2 = double(reshape(permute(enc(:,col1,sli(2:param.Rz/2:end),: ),[1 4 2 3]),param.N*param.ov*param.nCh,param.Rz*param.Ry/2));
%                         enc_res = [enc_res_1(:,1:param.Ry) enc_res_2(:,1:param.Ry) enc_res_1(:,param.Ry+1:end) enc_res_2(:,param.Ry+1:end)];

                for mm=1:param.Rz/param.caipi_del
                    aa = ii+shift_amount(mm):param.N/param.Ry:param.N+(param.N/param.Rz);
                    aa = aa(1:param.Ry);
                    col(mm,:) = aa;
                    if max(col(:)) > param.N
                        [val,idx] = max(col(:));
                        col(idx) = val-param.N;
                    end
                end
                col = repmat(col,param.caipi_del,1);
                for mm=1:param.Rz
                    enc_res(:,mm:param.Rz:end) = double(reshape(permute(enc(:,col(mm,:),sli(mm,:),:),[1 4 2 3]),param.N*param.ov*param.nCh,param.Rz));
                end
                for kk=1:param.N*param.ov
                        enc_res1 = enc_res(kk:param.N*param.ov:end,:);
                        EE = enc_res1'*enc_res1;
                        % Approach 1 Theoretical
                        EE_inv = pinv((EE));
                        EE_diag = diag(EE);
                        EE_inv_diag = (diag(EE_inv));
                        g_f = sqrt(EE_diag.*EE_inv_diag);
                        g_f = reshape(g_f,param.Ry,param.Rz);
                            for mm=1:param.Rz
                                g_img_t(kk,col(mm,:),sli(mm,:))=g_f(mm,:);
                            end
                end
                col = zeros(param.Ry/param.caipi_del,param.Rz);
        end
    end
       
% CAIPI unshifting
%         for ii=1:param.slices/param.Rz
%             g_img_t(:,:,ii:shift:param.slices,:) = circshift(g_img_t(:,:,ii:shift:param.slices,:),-param.N/2,2);
%         end
        
        %disp('Theoretical method: ')
        %toc

        msk=i_t(:,:,:,1);
        msk = padarray(msk,((param.N*param.ov)-param.N)/2,'both');
        msk(msk~=0) = 1;
        g_img_t = g_img_t.*msk;
        g_img_t = g_img_t((((param.N*param.ov)-param.N)/2)+1:(((param.N*param.ov)-param.N)/2)+param.N,:,:);
        g_av_t = abs(mean(g_img_t(g_img_t~=0)));
        g_max_t = abs(max(g_img_t,[],'all'));