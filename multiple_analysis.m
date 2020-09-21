% This is a bit of a messy script, but does the following:
%  1) find the wave-caipi parameters that are ok with the Slew Rate
%      and find values with similar RG, "params_use.params" has 
%      the similar configurations of the insert vs conventional systm
%  2) 
%     -   Loads the psf configurations for Equal/Different Gradient (ge/gd)
%         Equal/Different number of sines (se/de)
%     -  Plots the psf and Power spectrum for the selected configuration

%% 1) Finding acceptable WC-parameters
% Select parameters for the desired RG
des_rg = 1500;
% Manually tune this parameters to get values with similar Range of spread
tres_y = 10;
tres_z = 2;

N = 160;
p_banw = [50:50:600];
t_read = 1./p_banw;
G_sine = [1:2:40]./1000;
N_sine = [1:2:40];
SR_max = 200;
[T_R,G_S,N_S] = ndgrid(t_read,G_sine,N_sine);

SR_sine = 2.*pi.*(1./(T_R./N_S)).*G_S;
sinsy = N_S(SR_sine <=SR_max);
sinsz = sinsy;
gy = G_S(SR_sine <=SR_max);
gz = gy;
p_bw = 1./T_R(SR_sine <=SR_max);
sr = SR_sine(SR_sine <= SR_max);

datay = [gy gz sinsy sinsz p_bw];

% Calculate Range of spread
gamma=42.58e6;
t_r = 1./datay(:,5);
fy = datay(:,3)./t_r;
fz = datay(:,4)./t_r;
omgy = 2*pi.*fy;
omgz = 2*pi.*fz;
y1 = (-N/2)+1:N/2;
y1 = y1-mean(y1); y1 = y1.*1e-03;
z1=(-N/2)+1:N/2;
z1 = z1-mean(z1); z1 = z1.*1e-03;
t1=0:t_r/N:t_r-(t_r/N);
insfr = -gamma.*(datay(:,1).*cos(omgy.*t1).*y1+datay(:,2).*sin(omgz*t1).*z1);
fmax=max(insfr,[],2);
% rg=(2.*(0.001./param.p_bw).*fmax)+param.N;     % Range of spread
rg=(2./datay(:,5).*fmax)+N;     % Range of spread
datay = [datay rg];

% Calculating the time to get sines..
raise_t = gy./sr;
time = 4.*raise_t.*sinsy;
datay = [datay time];

figure; 
hold on
for idx = 1:length(p_banw)
    bw = p_banw(idx);
    aa = find(datay(:,5)==bw);
    aa = datay(aa,:);
    rg = aa(:,6);
    t = aa(:,7);
%     plot(bw,rg,'.','Color','b');
    subplot(1,2,1); plot(bw,rg,'.','Color','b');
    hold on
end
subplot(1,2,1);ylim([0 16000]); grid on

% Insert coil
p_bw = [50:50:600];
t_read = 1./p_bw;
G_sine = [6:4:40]./1000;
G_sine_is = [10:10:200]./1000;
N_sine = [5:12];
N_sine_is = [5:30];
SR_max = 200;
SR_max_is = 1300;
[T_R,G_S,N_S,G_is_S,N_is_S] = ndgrid(t_read,G_sine,N_sine,G_sine_is,N_sine_is);
SR_sine = 2.*pi.*(1./(T_R./N_S)).*G_S;
SR_sine_is = 2.*pi.*(1./(T_R./N_is_S)).*G_is_S;

index_valid_stuff = (SR_sine <=SR_max)&(SR_sine_is <=SR_max_is);
N_sine2use = N_S(index_valid_stuff);
N_is_sine2use = N_is_S(index_valid_stuff);
G_is_2use = G_is_S(index_valid_stuff);
G_2use = G_S(index_valid_stuff);
SR_sine_is_2use = SR_sine_is(index_valid_stuff);

P_bw_2use = 1./T_R(index_valid_stuff);

% Calculate Range of spread
gamma=42.58e6;
t_r = 1./P_bw_2use;
fy = N_sine2use./t_r;
fz = N_is_sine2use./t_r;
omgy = 2*pi.*fy;
omgz = 2*pi.*fz;
y1 = (-N/2)+1:N/2;
y1 = y1-mean(y1); y1 = y1.*1e-03;
z1=(-N/2)+1:N/2;
z1 = z1-mean(z1); z1 = z1.*1e-03;
t1=0:t_r/N:t_r-(t_r/N);
insfr = -gamma.*(G_2use.*cos(omgy.*t1).*y1+G_is_2use.*sin(omgz*t1).*z1);
fmax=max(insfr,[],2);
% rg=(2.*(0.001./param.p_bw).*fmax)+param.N;     % Range of spread
rg1=(2./P_bw_2use.*fmax)+N;     % Range of spread
sr_y = 2.*pi.*(1./(t_r./N_sine2use)).*G_2use;
sr_z = 2.*pi.*(1./(t_r./N_is_sine2use)).*G_is_2use;
rt_y = G_2use./sr_y;
rt_z = G_is_2use./sr_z;
ty = 4.*rt_y.*N_sine2use;
tz = 4.*rt_z.*N_is_sine2use;

data_ins = [G_2use G_is_2use N_sine2use N_is_sine2use P_bw_2use rg1 ty];

smple = [G_2use G_is_2use N_sine2use N_is_sine2use P_bw_2use rg1 ty SR_sine_is_2use];

% figure; plot(P_bw_2use,rg1,'o','Color','r'); grid on

subplot(1,2,2);  plot(P_bw_2use,rg1,'.','Color','r'); grid on

figure; plot(datay(:,6),datay(:,7),'o')
hold on
plot(data_ins(:,6),data_ins(:,7),'.')
xlim([0 4000]); grid on
 xlim([0 4000])
grid on

% Find the parameters with similar RG Insert vs Conventional Grad
diff_y = datay(:,6)-des_rg;
diff_z = data_ins(:,6)-des_rg;

y_use = abs(diff_y)<=tres_y;
z_use = abs(diff_z)<=tres_z;

aa_zeros = zeros([size(y_use(y_use==1),1) 1]);
bb_ones = ones([size(z_use(z_use==1),1) 1]);

aa = [datay(y_use,:) aa_zeros];
bb = [data_ins(z_use,:) bb_ones];

diff =  bb(:,5)-min(aa(:,5));
bb_i = diff>0;
bb = bb(bb_i,:);

params_use = [];
params_use.rg = des_rg;
params_use.params = [aa;bb];

tit = sprintf('./Data/rg_analysis/params_%i_%i.mat',des_rg,N);
% save(tit,'params_use')

clearvars -except aa bb params_use

%% 2) PSF analysis
% select parameter to plot psf and power
% res_gd_sd.param  -> gradient different, sines different
% res_gd_se.param  -> gradient different, sines equal
% res_ge_sd.param  -> gradient equal, sines different
% res_ge_se.param  -> gradient equal, sines equal
param = res_gd_sd.param;

load('./Data/psf_analysis/res_gd_sd.mat')
res_gd_sd  = res;
load('./Data/psf_analysis/res_gd_se.mat')
res_gd_se = res;
load('./Data/psf_analysis/res_ge_sd.mat')
res_ge_sd = res;
load('./Data/psf_analysis/res_ge_se.mat')
res_ge_se = res;

res_ge_se.pixels = sum(abs(res_ge_se.psf_yz(:))>1);
res_gd_se.pixels = sum(abs(res_gd_se.psf_yz(:))>1);
res_ge_sd.pixels = sum(abs(res_ge_sd.psf_yz(:))>1);
res_gd_sd.pixels = sum(abs(res_gd_sd.psf_yz(:))>1);

% Analyzing the power spectrum
psf_f=FFT_1D_Edwin(res_gd_sd.psf_yz,'image',1);
psf_f1 = psf_f(:,:,1);
t_r = 1/param.p_bw;
n = param.ov*param.Nx;
fq = (n)./t_r;
x=psf_f1(:,1);
x = FFT_1D_Edwin(x,'kspace',2);
x=x(:);
f=(-n/2:n/2-1)*(fq/n);
power=abs(x).^2/n;

%%% Plotting Power Spectrum, and psf
figure;subplot(2,1,1); imagesc(abs(psf_f1.')); title('Psf')
subplot(2,1,2); plot(f,power);  title('Power spectrum'); xlabel('Frequency Hz'); xlim([min(f) max(f)]);
hold on;
aa=0:0.001:max(power(:))+(0.25*max(power(:)));
fmax_four = max(f(power>(0.05*max(power(:)))));
bb=repmat(fmax_four,[1 size(aa,2)]);
plot(bb,aa,'LineWidth',1.5,'DisplayName','Max Fourier Frequency');
aa=0:0.001:max(power(:))+(0.25*max(power(:)));
bb=repmat(param.if_max,[1 size(aa,2)]);
plot(bb,aa,'LineStyle','--','LineWidth',1.5,'DisplayName','Max instant Frequency');legend;
colormap 'gray'
