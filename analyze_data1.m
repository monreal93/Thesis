% This is a messy script...
% Loads the data files from simulations and plots them, results
% are in folder "/Results/Scanner_CoilSens"

clear all ; clc

% Select parameters to load file
N=90;
slices=90;
Ry=3;
Rz=3;

tit = sprintf('./Results/Scanner_CoilSens/Data_%i_%i_%i_%i.mat',N,slices,Ry,Rz);

% tit = sprintf('/home/amonreal/Documents/Thesis/GitHub/Thesis/Results/Scanner_CoilSens/Data_160_80_14_14_70_4_4.mat');

data =  load(tit);
data = data.data;    

 %% Range of spread vs G-av
 % Generate index (weigthed sum of G-av, G-max and RMSE)
n_gav = (data(:,6)-min(data(:,6)))./(max(data(:,6))-min(data(:,6)));            % Normalized G-av
n_gmax = (data(:,7)-min(data(:,7)))./(max(data(:,7))-min(data(:,7)));       % Normalized G-max
n_rmse = (data(:,8)-min(data(:,8)))./(max(data(:,8))-min(data(:,8)));       % Normalized RMSE 
idx = (n_gav*0.4)+(n_gmax*0.2)+(n_rmse*0.4);
data = [data idx];

% Calculate Range of spread
gamma=42.58e6;
t_r = 1./data(:,5);
fy = data(:,3)./t_r;
fz = data(:,4)./t_r;
omgy = 2*pi.*fy;
omgz = 2*pi.*fz;
y1 = (-N/2)+1:N/2;
y1 = y1-mean(y1); y1 = y1.*data(:,10);
z1=(-N/2)+1:N/2;
z1 = z1-mean(z1); z1 = z1.*data(:,11);
t1=0:t_r/N:t_r-(t_r/N);
insfr = -gamma.*(data(:,1).*cos(omgy.*t1).*y1+data(:,2).*sin(omgz*t1).*z1);
fmax=max(insfr,[],2);
% rg=(2.*(0.001./param.p_bw).*fmax)+param.N;     % Range of spread
rg=(2./data(:,5).*fmax)+N;     % Range of spread
n_rg = (rg-min(rg))./(max(rg)-min(rg));   % Normalized Range of spread

% Rearenging Data matrix
data = [data rg n_rg];
data = sortrows(data,14);
idx = data(:,13); rg = data(:,14); n_rg = data(:,15);
g_av = data(:,6); g_max = data(:,7); rmse = data(:,8);

% % Plotting RG vs G-av/G-max/RMSE in separate plot
% f1= figure('Position',[-262,1151,1763,450]);
% tit = sprintf('Range of spread vs RMSE %ix%i',N,slices);
% subplot(1,3,1); scatter(rg,rmse); title(tit); xlabel('Spread out Range'); ylabel('RMSE');
% tit = sprintf('Range of spread vs G-av %ix%i',N,slices); 
% subplot(1,3,2); scatter(rg,g_av); title(tit); xlabel('Spread out Range'); ylabel('G-av');
% tit = sprintf('Range of spread vs G-max %ix%i',N,slices); 
% subplot(1,3,3); scatter(rg,g_max); title(tit); xlabel('Spread out Range'); ylabel('G-max');
% tit = sprintf('/home/amonreal/Documents/Thesis/Report/Figures/S10_results_%i_%i',N,slices);
% saveas(f1,tit,'epsc');


















