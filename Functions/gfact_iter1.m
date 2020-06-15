function [g_mean]=gfact_iter1(psf_yz,CoilSensitivity,param)

psf = repmat(psf_yz,[1 1 1 param.nCh]);
CoilSensitivity = padarray(CoilSensitivity,((param.N*param.ov)-param.N)/2,'both');
aa = (FFT_1D_Edwin(CoilSensitivity,'kspace',1));
enc = psf.*aa;
enc = FFT_1D_Edwin(enc,'image',1);

col = zeros(param.Ry/param.caipi_del,param.Rz);
enc_res = zeros(param.N*param.ov*param.nCh,param.Ry*param.Rz);

for mm=1:param.Rz-1
    shift_amount(mm+1) = mm;
end
shift_amount = shift_amount .* ((param.N/param.Rz)/(param.Rz/param.caipi_del));
shift_amount = floor(shift_amount);

% center slice with collapsed slices
sli = param.slices/2:param.slices/param.Rz:param.slices+(param.slices/param.Rz);
if max(sli(:)) > param.slices
    [val,idx] = max(sli(:));
    sli(idx) = val-param.slices;
end
sli=unique(sli);
sli = reshape(sli,param.Rz,[]);

% center colum with collapsed colums
for mm=1:param.Rz/param.caipi_del
    aa = (param.N/2)+shift_amount(mm):param.N/param.Ry:param.N+((param.N/param.Rz)*2);
    aa = aa(1:param.Ry);
    col(mm,:) = aa;
    if max(col(:)) > param.N
        idx = find(col>param.N);
        col(idx) = col(idx)-param.N;
    end
end
col = repmat(col,param.caipi_del,1);
                
for mm=1:param.Rz
    enc_res(:,mm:param.Rz:end) = double(reshape(permute(enc(:,col(mm,:),sli(mm,:),:),[1 4 2 3]),param.N*param.ov*param.nCh,param.Rz));
end
enc_res = enc_res((param.N*param.ov)/2:param.N*param.ov:end,:);
 
s = CoilSensitivity(:,:,param.N/2,:);
s = double(reshape(permute(s,[1 4 2 3]),param.N*param.ov*param.nCh,param.N));
SS = s'*s;

EE = enc_res'*enc_res;

% % Iterative approach:
% d = zeros(param.Ry*param.Rz,1);
% ec = zeros(param.Ry*param.Rz,1);
% ec((param.Ry*param.Rz/2),1) = 1;
% [d,flag]=lsqr(EE,ec);
% g_c = sqrt((ec'*d).*SS(param.N/2,param.N/2));

% Theroetical approach:
EE_inv = pinv((EE));
EE_diag = diag(EE);
EE_inv_diag = (diag(EE_inv));
g_f = sqrt(EE_diag.*EE_inv_diag);
g_f = reshape(g_f,param.Ry,param.Rz);
g_c = real(g_f(1));

g_mean = real(1+(0.37*(g_c-1)));

end

