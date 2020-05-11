function [g_img_t,g_av_t,g_max_t]=gfact_teo(N,slices,Ry,Rz,nCh,i_t,psf_yz,CoilSensitivity,ov,caipi)

        tic
        
        psf = repmat(psf_yz,[1 1 1 nCh]);
        aa = (FFT_1D_Edwin(CoilSensitivity,'kspace',1));
        enc = psf.*aa;
        enc = FFT_1D_Edwin(enc,'image',1);
        % 2D CAIPI effect
        enc_test1 = enc;
%         enc_test1 = circshift(enc_test1,N/Rz/2,3);
        for ii=1:slices/Rz
            if Rz==2
                shift=slices-(slices/(Rz*Rz));
            elseif Rz==4
                shift=(slices/Rz)*2;
            elseif Rz==1
                shift=1;
            end
            enc_test1(:,:,ii:shift:slices,:) = circshift(enc_test1(:,:,ii:shift:slices,:),N/2,2);
        end
        enc = enc_test1;
        
% Only collapsed colums and slices, pixel by pixel
g_img_t = zeros(N,N,slices);
   for jj=1:slices/Rz 
        sli = jj:slices/Rz:slices;
        for ii=1:N/Ry
                        col1=ii+(N/Ry/2):N/Ry:N+(N/Ry/2);
                        col2=ii:N/Ry:N; 
                        if max(col1) > N
                            col1(col1>N) = (col1(col1>N))-N;
                        end
                    
                        enc_res_1 = double(reshape(permute(enc(:,col2,sli(1:Rz/2:end-1),: ),[1 4 2 3]),N*ov*nCh,Rz*Ry/2));
                        enc_res_2 = double(reshape(permute(enc(:,col1,sli(2:Rz/2:end),: ),[1 4 2 3]),N*ov*nCh,Rz*Ry/2));
                        enc_res = [enc_res_1(:,1:Ry) enc_res_2(:,1:Ry) enc_res_1(:,Ry+1:end) enc_res_2(:,Ry+1:end)];
                    
            for kk=1:N*ov
                    enc_res1 = enc_res(kk:N*ov:end,:);
                    EE = enc_res1'*enc_res1;
                    % Approach 1 Theoretical
                    EE_inv = pinv((EE));
                    EE_diag = diag(EE);
                    EE_inv_diag = (diag(EE_inv));
                    g_f = sqrt(EE_diag.*EE_inv_diag);
                    g_f = reshape(g_f,1,Ry,Rz);
                    g_img_t(kk,col2,sli(1:Rz/2:end-1))=g_f(:,:,1:2:end);
                    g_img_t(kk,col1,sli(2:Rz/2:end))=g_f(:,:,2:2:end);
 
            end
        end
    end
       
% CAIPI unshifting
        for ii=1:slices/Rz
            g_img_t(:,:,ii:shift:slices,:) = circshift(g_img_t(:,:,ii:shift:slices,:),-N/2,2);
        end
        
        disp('Theoretical method: ')
        toc

        msk=i_t;
        msk = padarray(msk,((N*ov)-N)/2,'both');
        msk(msk~=0) = 1;
        g_img_t = g_img_t.*msk;
        g_img_t = g_img_t((((N*ov)-N)/2)+1:(((N*ov)-N)/2)+N,:,:);
        g_av_t = abs(mean(g_img_t(g_img_t~=0)));
        g_max_t = abs(max(g_img_t,[],'all'));
%         as(g_img_t)