% Script to calculate the G-factor values using the theoretical method, it
% first generates the whole encoding matrix (PSF*CoilSensitivity) and then
% subsets this matrix to the values that would collapse
% Note: Still need to check how to properly generate G-factor map (taking
% into account the oversampling factor)

function [g_img_t,g_av_t,g_max_t]=gfact_teo1(psf_yz,CoilSensitivity,param)

        psf = repmat(psf_yz,[1 1 1 param.nCh]);
        msk = CoilSensitivity;
%         CoilSensitivity = padarray(CoilSensitivity,((param.N*param.ov)-param.N)/2,'both');
        CoilSensitivity = padarray(CoilSensitivity,(size(psf,1)-size(CoilSensitivity,1))/2,'both');

        aa = (FFT_1D_Edwin(CoilSensitivity,'kspace',1));
        enc = psf.*aa;
        enc = FFT_1D_Edwin(enc,'image',1);
        
        % App 1
        enc = enc((((param.Nx*param.ov)-param.Nx)/2)+1:(((param.Nx*param.ov)-param.Nx)/2)+param.Nx,:,:,:);
        
% Only collapsed colums and slices, pixel by pixel
g_img_t = zeros(param.Nx,param.Ny,param.slices);
col = zeros(floor(param.Rz/param.caipi_del),param.Ry);
enc_res = zeros(param.Nx*param.nCh,param.Ry*param.Rz);

shift_amount = 0;
for mm=1:param.Rz-1
    shift_amount(mm+1) = mm;
end
shift_amount = shift_amount .* ((param.Ny/param.Rz)/(param.Rz/param.caipi_del));
shift_amount = floor(shift_amount);

   for jj=1:param.slices/param.Rz 
        sli = jj:param.slices/param.Rz:param.slices;
        sli = reshape(sli,param.Rz,[]);
        for ii=1:param.Ny/param.Ry
                        
                for mm=1:param.Rz/param.caipi_del
                    aa = ii+shift_amount(mm):param.Ny/param.Ry:param.Ny+(param.Ny/param.Rz);
                    aa = aa(1:param.Ry);
                    col(mm,:) = aa;
                    if max(col(:)) > param.Ny
                        idx = find(col>param.Ny);
                        col(idx) = col(idx)-param.Ny;
                    end
                end
                col = repmat(col,param.caipi_del,1);
                for mm=1:param.Rz
                    enc_res(:,mm:param.Rz:end) = double(reshape(permute(enc(:,col(mm,:),sli(mm,:),:),[1 4 2 3]),param.Nx*param.nCh,param.Ry));
                end
                for kk=1:param.Nx%*param.ov % App 2
                    enc_res1 = enc_res(kk:param.Nx:end,:);  % App 1
%                         enc_res1 = enc_res(kk:param.Nx*param.ov:end,:);  % App 2
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
   
   % App 2
%         g_img_t = g_img_t((size(g_img_t,1)/2)-(param.Nx/2):(size(g_img_t,1)/2)+(param.Nx/2)-1,:,:);
        msk= rssq(msk,4);
        msk(msk~=0) = 1;
        g_img_t = g_img_t.*msk;
        g_av_t = mean(abs(g_img_t(g_img_t~=0)));
        g_max_t = abs(max(g_img_t,[],'all'));
        
        