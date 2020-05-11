function [g_img_t,g_av_t,g_max_t]=gfact_teo1(N,slices,Ry,Rz,nCh,i_t,psf_yz,CoilSensitivity,ov)

        tic
        
        psf = repmat(psf_yz,[1 1 1 nCh]);
        % aa = (FFT_1D_Edwin(coilSensNew,'kspace',1));
        aa = (FFT_1D_Edwin(CoilSensitivity,'kspace',1));
        enc = psf.*aa;
        enc = FFT_1D_Edwin(enc,'image',1);
        enc = reshape(permute(enc,[1 4 2 3]),[],N,N);
        cs = reshape(permute(CoilSensitivity,[1 4 2 3]),[],N,N);
        clear psf aa psf_yz
        
        
        % Only collapsed colums and slices, pixel by pixel
        g_img_t = zeros(N*ov,N,slices);
        e_rho = zeros(N*N,1);
        d = zeros(N*N,1);
         
        for ii=1:N
            enc1 = enc(:,:,ii);
            cs1 = cs(:,:,ii);
            enc1 = reshape(enc1,nCh,[]);
            cs1 = reshape(cs1,nCh,[]);
            EE = enc1'*enc1;
            SS = cs1'*cs1;
            e_rho(ii,1) = 1;
            d = lsqr(EE,e_rho);
            
            g_f = sqrt((e_rho'*d)*diag(SS));
            g_f = reshape(g_f,N*ov,N);
            g_img_t(:,:,ii) = g_f;
            
            e_rho(ii,1) = 0;
        end
       
        
        disp('Theoretical method: ')
        toc

        msk=i_t;
        msk(msk~=0) = 1;
        g_img_t = g_img_t.*msk;
        g_av_t = abs(mean(g_img_t(g_img_t~=0)));
        g_max_t = abs(max(g_img_t,[],'all'));
%         as(g_img_t)