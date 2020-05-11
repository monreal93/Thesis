function [g_mean]=gfact_iter1(N,slices,Ry,Rz,nCh,psf_yz,CoilSensitivity,ov)

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

if Rz ==2
    sli=[32 64];
    col1=16+(N/Ry/2):N/Ry:N+(N/Ry/2);
    col2=16:N/Ry:N; 
elseif Rz==4
    sli=[16 32 48 64];
    col1=8+(N/Ry/2):N/Ry:N+(N/Ry/2);
    col2=8:N/Ry:N; 
end

enc_res_1 = double(reshape(permute(enc(:,col2,sli(1:Rz/2:end-1),: ),[1 4 2 3]),N*ov*nCh,Rz*Ry/2));
enc_res_2 = double(reshape(permute(enc(:,col1,sli(2:Rz/2:end),: ),[1 4 2 3]),N*ov*nCh,Rz*Ry/2));
enc_res = [enc_res_1(:,1:Ry) enc_res_2(:,1:Ry) enc_res_1(:,Ry+1:end) enc_res_2(:,Ry+1:end)];
s = CoilSensitivity(:,:,N/2,:);
s = double(reshape(permute(s,[1 4 2 3]),N*ov*nCh,N));
SS = s'*s;

EE = enc_res'*enc_res;

d = zeros(Ry*Rz,1);
ec = zeros(Ry*Rz,1);
ec((Ry*Rz/2),1) = 1;
[d,flag]=lsqr(EE,ec);

g_c = sqrt((ec'*d).*SS(N/2,N/2));
g_mean = real(1+(0.37*(g_c-1)));

