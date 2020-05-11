function i_wc = imag_wc(N,slices,Ry,Rz,caipi,caipi_del,kspace_nufft,plt)

if caipi && Ry ~= 1 || Rz ~= 1    

kspace_new = single(zeros([size(kspace_nufft,1) size(kspace_nufft,2)*Ry size(kspace_nufft,3)*Rz size(kspace_nufft,4)]));
kz_i=1;
ky_i=1;

% This one work for R3, I think this is the right one
        for ii=1:Rz
            kspace_new(:,ky_i:Ry*Ry:end,kz_i:Rz:end,:) = kspace_nufft(:,ii:Ry:end,:,:);  % for R=2
            kz_i = kz_i+caipi_del;
            ky_i = ky_i+(Ry);
            if kz_i>Rz
                kz_i=1;
            end
        end

    % Cropping from start
    lwy = 1; upy=N/Ry;
    lwz = 1; upz=slices/Rz;

    test_img = FFT_3D_Edwin(kspace_new,'image' );
    test_img1 = test_img(:,lwy:upy,lwz:upz,:);
    i_wc = test_img1;
else
    i_wc = FFT_3D_Edwin(kspace_nufft,'image' );
end