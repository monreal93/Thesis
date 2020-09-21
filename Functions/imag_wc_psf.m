% Generate the wave-CAIPI image, using the definition: psf and CoilSens

function imag_wc_psf = imag_wc_psf(CoilSensitivity,psf_yz,i_t,param)

  psf = repmat(psf_yz,[1 1 1 param.nCh]);
%   CoilSensitivity = CoilSensitivity.*i_t;
%   CoilSensitivity = padarray(CoilSensitivity,((param.N*param.ov)-param.N)/2,'both');      
%   aa = (FFT_1D_Edwin(CoilSensitivity,'kspace',1));
%   enc = psf.*aa;
  
  i_t =  padarray(i_t,((param.Nx*param.ov)-param.Nx)/2,'both');
  aa = (FFT_1D_Edwin(i_t,'kspace',1));
  enc = psf.*aa;
  enc = FFT_1D_Edwin(enc,'image',1);

  enc_f = FFT_3D_Edwin(enc,'kspace');

% Creating mask for CAIPI undersampling
msk = zeros(size(enc_f));
kz_i=1;
ky_i=1;
        for ii=1:param.Rz
                msk(:,ky_i:param.Ry*param.Rz:end,kz_i:param.Rz:end,:) = 1;
                kz_i = kz_i+param.caipi_del;
                ky_i = ky_i+param.Ry;
            if kz_i>param.Rz
                kz_i=1;
            end
        end
 

% Adding noise to k-space data
enc_f = addGaussianNoise_cmplx(enc_f,20); % Original 20

% openUndersampling in CAIPI way
enc_f = enc_f.*msk;
imag_wc_psf = FFT_3D_Edwin(enc_f,'image');
% Cropping from start
lwy = 1; upy=param.Ny/param.Ry;
lwz = 1; upz=param.slices/param.Rz;
imag_wc_psf = imag_wc_psf(:,lwy:upy,lwz:upz,:);
end