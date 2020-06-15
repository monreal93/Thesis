function i_wc = imag_wc(kspace_nufft,param)

if param.caipi && param.Ry ~= 1 || param.Rz ~= 1    

kspace_new = single(zeros([size(kspace_nufft,1) size(kspace_nufft,2)*param.Ry size(kspace_nufft,3)*param.Rz size(kspace_nufft,4)]));
kz_i=1;
ky_i=1;

        for ii=1:param.Rz
            if param.Ry == 1
                kspace_new(:,ky_i:(param.Ry*param.Rz)+(param.Rz-param.caipi_del):end,kz_i:param.Rz:end,:) = kspace_nufft(:,ii:param.Ry+(param.Rz-param.caipi_del):end,:,:);
            else
                kspace_new(:,ky_i:param.Ry*param.Rz:end,kz_i:param.Rz:end,:) = kspace_nufft(:,ii:param.Ry:end,:,:);
            end
                kz_i = kz_i+param.caipi_del;
            ky_i = ky_i+param.Ry;
            if kz_i>param.Rz
                kz_i=1;
            end
        end

    % Cropping from start
    lwy = 1; upy=param.N/param.Ry;
    lwz = 1; upz=param.slices/param.Rz;

    test_img = FFT_3D_Edwin(kspace_new,'image' );
    test_img1 = test_img(:,lwy:upy,lwz:upz,:);
    i_wc = test_img1;
else
    i_wc = FFT_3D_Edwin(kspace_nufft,'image' );
end
