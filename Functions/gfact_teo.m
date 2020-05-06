function [g_img_t,g_av_t,g_max_t]=gfact_teo(N,slices,Ry,Rz,nCh,i_t,psf_yz,CoilSensitivity,ov,caipi)

        tic
        
        psf = repmat(psf_yz,[1 1 1 nCh]);
        % aa = (FFT_1D_Edwin(coilSensNew,'kspace',1));
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
%         enc_test1 = circshift(enc_test1,-N/Rz/2,3);
        enc = enc_test1;
 
        %%
%         % Undersampling encoding time domain, suming collapsed ones
%         und_enc = zeros(N*ov,N/Ry,N/Rz,nCh);
%        
% %         for ii=1:N/Ry
% %             for jj=1:N/Rz
% %                 enc_ff = enc(:,jj:N/Ry:end,ii:N/Rz:end,:);
% %                 enc_ff = sum(sum(enc_ff,2),3);
% %                 und_enc(:,ii,jj,:) = enc_ff;
% %             end
% %         end
%         
%    und_enc = FFT_3D_Edwin(und_enc,'kspace');
%       % Getting CAIPI pattern, adding zeros
%          % Define the start position for the shifted row in CAIPI 
%     if Rz==4
%       mm=3;
%     else
%       mm=2;
%     end
%  kk=1;
%  enc_und_full = zeros(size(enc*Ry));
%     for iter = 1:N
%             if mod(iter,2) == 0
%                 enc_und_full(:,kk,mm:Rz:end,:) = und_enc(:,iter,:,:);
%                 kk = kk+Ry;
%             else
%                  enc_und_full(:,kk,1:Rz:end,:) = und_enc(:,iter,:,:);
%                  kk = kk+Ry;
%             end 
%     end
%       enc_und_full = FFT_3D_Edwin(enc_und_full,'image');  
%       enc = enc_und_full;  
        
 %%     
%         enc2 = zeros(size(enc));
%         
%         if Rz == 1
%             shift_z = 0;
%         elseif Rz ==2
%             shift_z = [-1,0] .* N/Rz/ 2;
%         elseif Rz ==3
%             shift_z = [-1,0,1] .* N/Rz/ 2;
%         elseif Rz == 4
%             shift_z = [-2,-1,0,1] .* N/Rz/ 2;
%         end
%         if Ry == 1
%             shift_y = 0;
%         elseif Ry ==2
%             shift_y = [-1,0] .* N/Ry/2;
%         elseif Ry ==3
%             shift_y = [-1,0,1] .* N/Ry/2;
%         elseif Ry == 4
%             shift_y = [-2,-1,0,1] .* N/Ry/2;
%         end
%         
%         for nz = 1:Ry
% %             enc2(:,:,nz,:) = enc2(:,:,nz,:) +circshift(enc(:,:,nz,:), [0,shift_y(nz),0,0]);
% %             enc2(:,:,nz:Ry:end,:) = enc2(:,:,nz:Ry:end,:) +circshift(enc(:,:,nz:Ry:end,:), [0,shift_y(nz),0,0]);
%                     enc2 = enc2+circshift(enc,[0 shift_y(nz) 0 0]);
%         end
%         
% %         enc2 = circshift(enc2,[0 0 shift_z(1) 0]);
% %         for nz = 1:Rz
% %             enc2 = circshift(enc, [0,0,shift_z(nz),0]);
% %         end
%         % To shift y
%         enc2 = circshift(enc,[0 shift_y(1) shift_z(2) 0]);
%         enc3 = circshift(enc,[0 -shift_y(1) shift_z(2) 0]);
%         % To shift z
%         enc4 = circshift(enc,[0 shift_y(1) shift_z(1) 0]);
%         enc5 = circshift(enc,[0 shift_y(1) -shift_z(1) 0]);
%         
%         final_enc = enc + enc2+enc3+enc4+enc5;
%         enc = final_enc;

%%
        % Only collapsed colums and slices, pixel by pixel
        g_img_t = zeros(N,N,slices);
        for jj=1:slices/Rz  %2D
        sli = jj:slices/Rz:slices;
            for ii=1:N/Ry
%                     col = ii:N/Ry:N;
                    % CAIPI
                    %%%%% NEED TO ADD THIS IN A LOOP, to work with R>2
                    if Ry == 2
                        col1=ii+(N/Ry/2):N/Ry:N+(N/Ry/Rz);
                        col2=ii:N/Ry:N; 
                        if max(col1) > N
                            col1(col1>N) = (col1(col1>N))-N;
                        end
                    end
                    
                    if Ry == 4
                        col1=ii+(N/Ry/2):N/Ry:N+(N/Ry/Rz)+5;
    %                     col1=ii:N/Ry:N;
                        col2=ii:N/Ry:N; 
                        if max(col1) > N
                            col1(col1>N) = (col1(col1>N))-N;
                        end
                    end
                    
                    if Rz==1
                        enc_res = double(reshape(permute(enc(:,col2,sli,: ),[1 4 2 3]),N*ov*nCh,Rz*Ry));
                    else
                        if caipi~= 1
                            col1=col2;
                        end
                        enc_res_1 = double(reshape(permute(enc(:,col2,sli(1:Rz/2:end-1),: ),[1 4 2 3]),N*ov*nCh,Rz*Ry/2));
                        enc_res_2 = double(reshape(permute(enc(:,col1,sli(2:Rz/2:end),: ),[1 4 2 3]),N*ov*nCh,Rz*Ry/2));
                        enc_res = [enc_res_1(:,1:Ry) enc_res_2(:,1:Ry) enc_res_1(:,Ry+1:end) enc_res_2(:,Ry+1:end)];
                    end
                    
                    %No-CAIPI
%                     enc_res = double(reshape(permute(enc(:,col,sli,: ),[1 4 2 3]),N*ov*nCh,Ry*Rz));
% %                 enc_res = enc(:,col,sli,:);
% %                 enc_res = squeeze(sum(sum(enc_res,2),3));
% %                 enc_res = enc_res(:);
% %                 enc_res_test = sum(enc_res,2);   % New
                for kk=1:N*ov
                    enc_res1 = enc_res(kk:N*ov:end,:);
%                     enc_res1 = sum(enc_res1,1);
                    EE = enc_res1'*enc_res1;
%                     EE = EE./(nCh-1);
                    % Approach 1 Theoretical
                    EE_inv = pinv((EE));
                    EE_diag = diag(EE);
                    EE_inv_diag = (diag(EE_inv));
                    g_f = sqrt(EE_diag.*EE_inv_diag);
                    % CAIPI
                    g_f = reshape(g_f,1,Ry,Rz);
                    
                    if Rz==1
                        g_img_t(kk,col2,sli)=g_f;
                    else
                        g_img_t(kk,col2,sli(1:Rz/2:end-1))=g_f(:,:,1:2:end);
                        g_img_t(kk,col1,sli(2:Rz/2:end))=g_f(:,:,2:2:end);
                    end

                    % No-CAIPI
%                     g_f = reshape(g_f,1,Ry,Rz);
%                     g_img_t(kk,col,sli)=g_f(:,:,:);
                end
            end
        end
       
% CAIPI unshifting
%         g_img_t = circshift(g_img_t,N/Rz/2,3);
        for ii=1:slices/Rz
            g_img_t(:,:,ii:shift:slices,:) = circshift(g_img_t(:,:,ii:shift:slices,:),-N/2,2);
        end
%         g_img_t = circshift(g_img_t,-N/Rz/2,3);
        
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