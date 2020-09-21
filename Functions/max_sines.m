% Just calculating the Max number of sines that can be used with the given
% configuration, if selected param is larger, the max one would be used

function param = max_sines(param)
    t_r = 1./param.p_bw;
    gz_rt = 1/(param.gz_sr/param.gz);           % Z Gradient raise time
    gy_rt = 1/(param.gy_sr/param.gy);          % Y Gradient raise time
    
    % Trying another approach , SR not fixed....
%     sr_y = 2.*pi.*(1./(t_r./param.sinsy)).*param.gy;
%     sr_z = 2.*pi.*(1./(t_r./param.sinsz)).*param.gz;
%     gz_rt = 1/(sr_z/param.gz);           % Z Gradient raise time
%     gy_rt = 1/(sr_y/param.gy);          % Y Gradient raise time

    max_sinsy = floor(t_r./(4*gy_rt));
    max_sinsz = floor(t_r./(4*gz_rt));
    if param.sinsy > max_sinsy
        param.sinsy = max_sinsy;
    end
    if param.sinsz > max_sinsz
        param.sinsz = max_sinsz;
    end
end