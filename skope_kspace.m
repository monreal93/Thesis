% clearvars;

%% Sources
addpath './/bin'
addpath(genpath('/home/amonreal/Documents/Thesis/Scanner_programming/Field Camera'))

%% Example data from a prototype insert-gradient system
dataFolder = '/home/amonreal/Documents/Thesis/GitHub/Thesis/Data/wc_scan_230620/WaveCaipiPSF';
scanId_y = 4;
scanId_z = 3;

%% read Y data
scan_y = AqSysData(dataFolder,scanId_y);
raw_y = scan_y.getData('raw', [], [], [], 1);
kcoco_y =         scan_y.getData('kcoco', [], [], [], 1);
kspha_y =         scan_y.getData('kspha', [], [], [], 1);

%% read Z data
scan_z = AqSysData(dataFolder,scanId_z);
raw_z = scan_z.getData('raw', [], [], [], 1);
kcoco_z =         scan_z.getData('kcoco', [], [], [], 1);
kspha_z =         scan_z.getData('kspha', [], [], [], 1);

%% Creating cartesian grid
param.N = 224;
param.slices = 60;
param.ov = 6;
param.p_s = [1 1 2]*0.001;
param.i_fov = [param.N param.N param.slices].*param.p_s;
% Creating cartesian coordinates
yc = param.i_fov(2)/2*-1:param.p_s(2):(param.i_fov(2)/2)-param.p_s(2);
yc = yc - mean(yc);
zc_2mm = param.i_fov(3)/2*-1:param.p_s(3):(param.i_fov(3)/2)-param.p_s(3);
zc_2mm = zc_2mm - mean(zc_2mm);
zc_1mm = yc;

%% Generating Z k-space
    p_bw = 70;
    t_r = 1./p_bw;
    dt = 9.300000034272795e-06;
    ds_w = 13.017600059509277e-3;
    dt_acq = t_r./224*1.37;
    fq = 1/ds_w/7;
    samples_shift = ((ds_w/7)/2)/dt;
    K_ramp = 42.577e6*0.5*0.006*(0.006/200)*2*pi;
    rmp = yc.*K_ramp;

    k_z = kspha_z(:,4);
    k_zz = k_z.*zc_1mm;

    T_fieldmeasurement = 0:dt:(dt.*(length(kspha_z)-1));
    t_acq = 0:dt_acq:dt_acq*(1344-1);
    start_aq = 3700*dt;

    Wave_z = interp1(T_fieldmeasurement,k_zz,t_acq+start_aq,'cubic',0);
%     Wave_z = Wave_z - mean(Wave_z);             % to make it sine like 
%     Wave_z = circshift(Wave_z,[50 0]);              % to make it sine like

    psf_z = exp(-1i*Wave_z);
    

%% Generating Y k-space 
    p_bw = 70;
    t_r = 1./p_bw;
    dt = 9.300000034272795e-06;
    ds_w = 13.017600059509277e-3;
    dt_acq = t_r./224;
    fq = 1/ds_w/7;
    samples_shift = ((ds_w/7)/2)/dt;
    K_ramp = 42.577e6*0.5*0.006*(0.006/200)*2*pi;
    rmp = yc.*K_ramp;
    
    k_y = kspha_y(:,2);
    k_yy = k_y.*yc;
    
    T_fieldmeasurement = 0:1e-6:(1e-6.*(length(kspha_y)-1));
    t_acq = 0:dt:dt*(1344-1);
    start_aq = 3700*1e-6; %3700
    
    Wave_y = interp1(T_fieldmeasurement,k_yy,t_acq+start_aq,'cubic',0);
%     Wave_y = Wave_y - mean(Wave_y);             % to make it sine like 
%     Wave_y = circshift(Wave_y,[50 0]);              % to make it sine like
    
    psf_y = exp(-1i*Wave_y);
    
    
%     [psf_y, psf_z] = deal(psf_z,psf_y);
    save('Data/psf_z.mat','psf_z');
    save('Data/psf_y.mat','psf_y');
