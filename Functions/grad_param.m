function [x,y,z,x_calc,y_calc,z_calc,del_ky,del_kz,t,t_prime] = grad_param(gy,gz,sinsy,sinsz,p_bw,p_s,k_fov,i_fov,Ry,Rz,x_points,gamma)
    
t_r = 1/p_bw;
delta_t = 1e-7;
t_prime = 0:delta_t:t_r;
t = 0:t_r/x_points:t_r;
fy = sinsy/t_r;
fz = sinsz/t_r;
omgy = 2*pi*fy;
omgz = 2*pi*fz;
x_calc_points = length(t_prime);

%For helix....
x_calc = ((-k_fov(1)/2)+(k_fov(1)/x_calc_points)):(k_fov(1))/x_calc_points:(k_fov(1)/2);
x_calc = x_calc-mean(x_calc);
y_calc = gamma*cumsum(gy*cos(omgy*t_prime)*delta_t);
%y_calc = y_calc - mean(y_calc);
z_calc = gamma*cumsum(gz*sin(omgz*t_prime)*delta_t);
%z_calc = z_calc - mean(z_calc);
range_y = (k_fov(1)/2)-max(y_calc);
range_z = (k_fov(1)/2)-max(z_calc);
x = interp1(t_prime,x_calc,t,'linear','extrap');
y = interp1(t_prime,y_calc,t,'linear','extrap')+range_y;
z =  interp1(t_prime,z_calc,t,'linear','extrap')+range_z;

% Defining spacing between helix in y and z directions
del_ky = (1/(i_fov(2)*p_s(2))*Ry);
del_kz = (1/(i_fov(3)*p_s(3))*Rz);


