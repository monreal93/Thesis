% Calculates the RMSE, it uses as a reference the selected image "r_img",
% either phantom or a path of the reference volume

function [rmse] = norm_mse(r_img,i_wc_recon,CoilSensitivity,param)

if r_img == "phantom"
    r_img=phantom3d(param.Nx);
    
    if param.Nx ~= param.slices
        r_img = r_img(:,:,((param.Nx-param.slices)/2)+1:((param.Nx-param.slices)/2)+param.slices);
    end
    
else
    r_img = load(r_img);
    cells = struct2cell(r_img);
    r_img = cells{1};
   
end

msk = CoilSensitivity(:,:,:,1);
msk(msk~=0) = 1;

r_img = padarray(r_img,[(size(msk,1)-size(r_img,1))/2 ...
    (size(msk,2)-size(r_img,2))/2 (size(msk,3)-size(r_img,3))/2],'both');

r_img = double(r_img.*msk);
r_img = abs(r_img);
i_wc_recon = abs(i_wc_recon);

% Normalizing 
% i_wc_recon = (i_wc_recon-min(i_wc_recon(i_wc_recon~=0)))/(max(i_wc_recon(i_wc_recon~=0))-min(i_wc_recon(i_wc_recon~=0)));
% r_img = (r_img-min(r_img(r_img~=0)))/(max(r_img(r_img~=0))-min(r_img(r_img~=0)));

% Need to check how to properly scale it...
% rmse = sqrt(immse(i_wc_recon,r_img)./abs(mean(r_img(:)>eps)));

i_wc_recon(isnan(i_wc_recon))=0;
r_img(isnan(r_img))=0;
rse = ((abs((i_wc_recon)./median(i_wc_recon(r_img>eps))-r_img./median(r_img(r_img>eps))).^2))./(r_img./median(r_img>eps));
rse = rse.*(r_img>eps);

rmse =sqrt(mean(rse(r_img(:)>eps)));

% sqrt(abs((i_wc_recon)./median(i_wc_recon(r_img(:)>0))-r_img./median(r_img(r_img(:)>0)))).^2;

end


