function [g_img_i,g_av_i,g_max_i]=gfact_iter(N,slices,Ry,Rz,nCh,i_t,psf_yz,CoilSensitivity,ov)
        tic

        psf = repmat(psf_yz,[1 1 1 nCh]);
        % aa = (FFT_1D_Edwin(coilSensNew,'kspace',1));
        aa = (FFT_1D_Edwin(CoilSensitivity,'kspace',1));
        enc = psf.*aa;
        enc = FFT_1D_Edwin(enc,'image',1);
        % enc = CoilSensitivity;

%%%%%%%%%%%%%%%%%%% 2D CAIPI effect
        enc_test1 = enc;
%         enc_test1 = circshift(enc_test1,N/Rz/2,3);
        for ii=1:slices/Rz
            if Rz==2
                shift=slices-(slices/(Rz*Rz));
            elseif Rz==4
                shift=(slices/Rz)*2;
            end
            enc_test1(:,:,ii:shift:slices,:) = circshift(enc_test1(:,:,ii:shift:slices,:),N/2,2);
        end
%         enc_test1 = circshift(enc_test1,-N/Rz/2,3);
        enc = enc_test1;
%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Only collapsed colums and slices...
        g_img_i = zeros(N,N,slices);
        for jj=1:slices/Rz  %2D
        sli = jj:slices/Rz:slices;
            for ii=1:N/Ry
%                 col=ii:N/Ry:N;
%                 enc_res = single(reshape(permute(enc(:,col,sli,: ),[1 4 2 3]),N*ov*nCh,Ry*Rz));
                
%%%%%%%%% 2D CAIPI effect
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
                    
                    
                    enc_res_1 = double(reshape(permute(enc(:,col2,sli(1:Rz/2:end-1),: ),[1 4 2 3]),N*ov*nCh,Rz*Ry/2));
                    enc_res_2 = double(reshape(permute(enc(:,col1,sli(2:Rz/2:end),: ),[1 4 2 3]),N*ov*nCh,Rz*Ry/2));
                    enc_res = [enc_res_1(:,1:Ry) enc_res_2(:,1:Ry) enc_res_1(:,Ry+1:end) enc_res_2(:,Ry+1:end)];
                    
%%%%%%%%%%%
                
                
                for kk=1:N
                    enc_res1 = enc_res(kk:N*ov:end,:);
%                     enc_res1 = enc(kk,col,sli,:);
%                     enc_res1 = reshape(permute(enc_res1,[1 4 2 3]),nCh,[]);
                    EE = enc_res1'*enc_res1;
                    % Approach 2 iterative
                    g_f = zeros(1,size(EE,1));
                    for ll=1:size(EE,1)
                        e_rho = zeros(size(EE,1),1);
                        e_rho(ll) = 1;
                        [d,flag]=lsqr(EE,e_rho);
                        g_f(ll) = sqrt(e_rho'*d)*sqrt(e_rho'*EE*e_rho);
                    end
                        g_f = reshape(g_f,1,Ry,Rz);
%                         g_img_i(kk,col,sli)=g_f;
                    g_img_i(kk,col2,sli(1:Rz/2:end-1))=g_f(:,:,1:2:end);
                    g_img_i(kk,col1,sli(2:Rz/2:end))=g_f(:,:,2:2:end);
                end
            end
        end

%%%%%%%%% CAIPI unshifting
%         g_img_t = circshift(g_img_t,N/Rz/2,3);
        for ii=1:slices/Rz
            g_img_i(:,:,ii:shift:slices,:) = circshift(g_img_i(:,:,ii:shift:slices,:),-N/2,2);
        end
%         g_img_t = circshift(g_img_t,-N/Rz/2,3);

%%%%%%%%%%%%
        
        disp('Iterative method: ')
        toc

        msk=i_t;
        msk(msk~=0) = 1;
        g_img_i = g_img_i.*msk;
        g_av_i = mean(g_img_i(g_img_i~=0));
        g_max_i = max(g_img_i,[],'all');
        as(g_img_i)