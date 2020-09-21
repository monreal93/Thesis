% Calculates the oversampling factor, using the instantaneous frequency, it
% calcualtes the Range of Spread as well

function [ov,rg,if_max] = ov(param)
%% Calculating over sampling
gamma=42.58e6;
t_r = 1/param.p_bw;
fy = param.sinsy/t_r;
fz = param.sinsz/t_r;
omgy = 2*pi*fy;
omgz = 2*pi*fz;
y1 = (-param.Ny/2)+1:param.Ny/2;
y1 = y1-mean(y1); y1 = y1*param.p_s(2);
z1=(-param.Ny/2)+1:param.Ny/2;
z1 = z1-mean(z1); z1 = z1*param.p_s(3);
t1=0:t_r/param.Nx:t_r-(t_r/param.Nx);
insfr = -gamma*(param.gy*cos(omgy*t1).*y1+param.gz*sin(omgz*t1).*z1);
if_max=max(insfr,[],2);
rg=(2./param.p_bw.*if_max)+param.Nx;     % Range of spread

ov = ceil((rg/param.Nx));
ov_max = floor((t_r/(1e-6))./param.Nx);
if ov > ov_max
    ov = ov_max;
end

end