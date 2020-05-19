function CoilSensitivity = coil_sens(i_t,param)
[CoilSensitivity1] = GenerateCoilSensitivity4Sim(size(i_t),[param.Resolution param.Resolution param.Resolution],param.coil_radius,param.loop_radius,param.z_offset(1),param.nCh./4);
[CoilSensitivity2] = GenerateCoilSensitivity4Sim(size(i_t),[param.Resolution param.Resolution param.Resolution],param.coil_radius,param.loop_radius,param.z_offset(2),param.nCh./4);
[CoilSensitivity3] = GenerateCoilSensitivity4Sim(size(i_t),[param.Resolution param.Resolution param.Resolution],param.coil_radius,param.loop_radius,param.z_offset(3),param.nCh./4);
[CoilSensitivity4] = GenerateCoilSensitivity4Sim(size(i_t),[param.Resolution param.Resolution param.Resolution],param.coil_radius,param.loop_radius,param.z_offset(4),param.nCh./4);
CoilSensitivity = cat(4,CoilSensitivity1,CoilSensitivity2,CoilSensitivity3,CoilSensitivity4);

mask = ones(size(i_t));
mask(i_t == 0) = 0;
mask = repmat(mask,1,1,1,param.nCh);
CoilSensitivity = CoilSensitivity.*mask;
CoilSensitivity = CoilSensitivity./max(CoilSensitivity(:));

if param.plt == 1
    as(CoilSensitivity)
end