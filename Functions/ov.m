function ov = ov(param)
%% Calculating over sampling
gamma=42.58e6;
t_r = 1/param.p_bw;
fy = param.sinsy/t_r;
fz = param.sinsz/t_r;
omgy = 2*pi*fy;
omgz = 2*pi*fz;
y1 = (-param.N/2)+1:param.N/2;
y1 = y1-mean(y1);
z1=(-param.N/2)+1:param.N/2;
z1 = z1-mean(z1);
t1=0:t_r/param.N:t_r-(t_r/param.N);
insfr = -gamma*(param.gy*cos(omgy*t1).*y1+param.gz*sin(omgz*t1).*z1);
fmax=max(insfr);
rg=(2*(0.001/param.p_bw)*fmax)+param.N;
ov = ceil((rg/param.N));

end