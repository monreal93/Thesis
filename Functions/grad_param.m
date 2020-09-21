% This script calculates the corckscrew, saves the
% points in the space (taking into account the oversampling factor)
% and defines the disance between each corckscrew, 

function [x,y,z,x_calc,y_calc,z_calc,del_ky,del_kz,t,t_prime] = grad_param(param)
    
% t_r = 1/param.p_bw;
t_r = param.t_r;
delta_t = 1e-7;
t_prime = 0:delta_t:t_r;
t = 0:t_r/param.x_points:t_r;
fy = param.sinsy/t_r;
fz = param.sinsz/t_r;
omgy = 2*pi*fy;
omgz = 2*pi*fz;
x_calc_points = length(t_prime);

%For helix....
x_calc = ((-param.k_fov(1)/2)+(param.k_fov(1)/x_calc_points)):(param.k_fov(1))/x_calc_points:(param.k_fov(1)/2);
% x_calc = x_calc-mean(x_calc);
y_calc = param.gamma*cumsum(param.gy*cos(omgy*t_prime)*delta_t);
% y_calc = y_calc - mean(y_calc);
z_calc = param.gamma*cumsum(param.gz*sin(omgz*t_prime)*delta_t);
% z_calc = z_calc - mean(z_calc);
range_y = (param.k_fov(1)/2)-max(y_calc);
range_z = (param.k_fov(1)/2)-max(z_calc);
x = interp1(t_prime,x_calc,t,'linear','extrap');
y = interp1(t_prime,y_calc,t,'linear','extrap')+range_y;
z =  interp1(t_prime,z_calc,t,'linear','extrap')+range_z;

% Defining spacing between helix in y and z directions
del_ky = (1/(param.i_fov(2)*param.p_s(2))*param.Ry);
del_kz = (1/(param.i_fov(3)*param.p_s(3))*param.Rz);


