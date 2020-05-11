function CoilSensitivity = coil_sens(z_offset,coil_radius,loop_radius,Resolution,i_t,nCh,plt)

% [CoilSensitivity1] = GenerateCoilSensitivity4Sim(size(i_t),[Resolution Resolution Resolution],coil_radius,loop_radius,z_offset(1),nCh./2);
% [CoilSensitivity2] = GenerateCoilSensitivity4Sim(size(i_t),[Resolution Resolution Resolution],coil_radius,loop_radius,z_offset(2),nCh./2);
% CoilSensitivity = cat(4,CoilSensitivity1,CoilSensitivity2);

[CoilSensitivity1] = GenerateCoilSensitivity4Sim(size(i_t),[Resolution Resolution Resolution],coil_radius,loop_radius,z_offset(1),nCh./4);
[CoilSensitivity2] = GenerateCoilSensitivity4Sim(size(i_t),[Resolution Resolution Resolution],coil_radius,loop_radius,z_offset(2),nCh./4);
[CoilSensitivity3] = GenerateCoilSensitivity4Sim(size(i_t),[Resolution Resolution Resolution],coil_radius,loop_radius,z_offset(3),nCh./4);
[CoilSensitivity4] = GenerateCoilSensitivity4Sim(size(i_t),[Resolution Resolution Resolution],coil_radius,loop_radius,z_offset(4),nCh./4);
CoilSensitivity = cat(4,CoilSensitivity1,CoilSensitivity2,CoilSensitivity3,CoilSensitivity4);

mask = ones(size(i_t));
mask(i_t == 0) = 0;
mask = repmat(mask,1,1,1,nCh);
CoilSensitivity = CoilSensitivity.*mask;
CoilSensitivity = CoilSensitivity./max(CoilSensitivity(:));



if plt == 1
    as(CoilSensitivity)
end