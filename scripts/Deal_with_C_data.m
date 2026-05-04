clear; clc;
%%
nx1=48;
nx2=32;
nx3=64;

i = 2;
j = 2;
k = 1;

NO2 = 4;
ww = find_max_seq();

x1c = stitch_hdf5_array(i,j,k,NO2,'/x1c',0,0,0,1,1,0);
x2c = stitch_hdf5_array(i,j,k,NO2,'/x2c',0,0,0,1,1,0);
x3c = stitch_hdf5_array(i,j,k,NO2,'/x3c',0,0,0,1,1,0);
x1 = stitch_hdf5_array(i,j,k,NO2,'/x1',1,1,1,1,1,0);
x2 = stitch_hdf5_array(i,j,k,NO2,'/x2',1,1,1,1,1,0);
x3 = stitch_hdf5_array(i,j,k,NO2,'/x3',1,1,1,1,1,0);

x1ctr = stitch_hdf5_array(i,j,k,NO2,'/x1ctr',0,0,0,1,1,0);
x2ctr = stitch_hdf5_array(i,j,k,NO2,'/x2ctr',0,0,0,1,1,0);
x3ctr = stitch_hdf5_array(i,j,k,NO2,'/x3ctr',0,0,0,1,1,0);

dx1 = stitch_hdf5_array(i,j,k,NO2,'/dx1',0,1,1,1,1,0);
dx2 = stitch_hdf5_array(i,j,k,NO2,'/dx2',1,0,1,1,1,0);
dx3 = stitch_hdf5_array(i,j,k,NO2,'/dx3',1,1,0,1,1,0);

dxi = stitch_hdf5_array(i,j,k,NO2,'/geo',0,1,1,4,4,0);
dxj = stitch_hdf5_array(i,j,k,NO2,'/geo',1,0,1,5,5,0);
dxk = stitch_hdf5_array(i,j,k,NO2,'/geo',1,1,0,6,6,0);

dxi_ring = stitch_hdf5_array(i,j,k,NO2,'/geo',0,1,1,7,7,0);
dxj_ring = stitch_hdf5_array(i,j,k,NO2,'/geo',1,0,1,8,8,0);
dxk_ring = stitch_hdf5_array(i,j,k,NO2,'/geo',1,1,0,9,9,0);

face_idir = stitch_hdf5_array(i,j,k,NO2,'/geo',1,0,0,10,10,0);
face_jdir = stitch_hdf5_array(i,j,k,NO2,'/geo',0,1,0,11,11,0);
face_kdir = stitch_hdf5_array(i,j,k,NO2,'/geo',0,0,1,12,12,0);
vol = stitch_hdf5_array(i,j,k,NO2,'/geo',0,0,0,16,16,0);

B01 = stitch_hdf5_array(i,j,k,NO2,'/B0_i',1,0,0,1,1,0);
B02 = stitch_hdf5_array(i,j,k,NO2,'/B0_j',0,1,0,1,1,0);
B03 = stitch_hdf5_array(i,j,k,NO2,'/B0_k',0,0,1,1,1,0);

arrays_to_permute = {x1c, x2c, x3c, x1, x2, x3, x1ctr, x2ctr, x3ctr, ...
                    dx1, dx2, dx3, face_idir, face_jdir, face_kdir, vol, ...
                    B01, B02, B03};

array_names = {'x1c', 'x2c', 'x3c', 'x1', 'x2', 'x3', 'x1ctr', 'x2ctr', 'x3ctr', ...
              'dx1', 'dx2', 'dx3', 'face_idir', 'face_jdir', 'face_kdir', 'vol', ...
              'B01', 'B02', 'B03'};

for idx = 1:length(arrays_to_permute)
    current_array = arrays_to_permute{idx};
    numDims = ndims(current_array);
    newOrder = numDims:-1:1;
    eval([array_names{idx} ' = permute(current_array, newOrder);']);
end

[nx_total, ny_total, nz_total] = size(x1);

ic_act = NO2+1:NO2+nx1;
jc_act = NO2+1:NO2+nx2;
kc_act = NO2+1:NO2+nx3;

if_act = NO2+1:NO2+nx1+1;
jf_act = NO2+1:NO2+nx2+1;
kf_act = NO2+1:NO2+nx3+1;

x1ctr_act = x1ctr;
x2ctr_act = x2ctr;
x3ctr_act = x3ctr;

face_idir_act = face_idir;
face_jdir_act = face_jdir;
face_kdir_act = face_kdir;
vol_act = vol;

B01_act = B01;
B02_act = B02;
B03_act = B03;

% Change based on your conditions
normalization_params = [
    1.0;    % x_Norm
    1.0;    % u_Norm
    1.0;    % Time_Norm
    1.0;    % rho_Norm
    1.0;    % p_Norm
    1.0;    % B_Norm
];

for time_step = 0:ww
    fprintf('Processing time step %d of %d...\n', time_step, ww);
    
    try
        gas_mpi = stitch_hdf5_array(i,j,k,NO2,'/gas',0,0,0,1,5,time_step);
        mag1_mpi = stitch_hdf5_array(i,j,k,NO2,'/gem',1,0,0,1,1,time_step);
        mag2_mpi = stitch_hdf5_array(i,j,k,NO2,'/gem',0,1,0,2,2,time_step);
        mag3_mpi = stitch_hdf5_array(i,j,k,NO2,'/gem',0,0,1,3,3,time_step);
        
        filename = sprintf('../build//mhd_%02d-%02d-%02d_%06d.h5', 0, 0, 0, time_step);
        t = h5read(filename, '/time_sim');
        
        gas_mpi = permute(gas_mpi, ndims(gas_mpi):-1:1);
        mag1_mpi = permute(mag1_mpi, ndims(mag1_mpi):-1:1);
        mag2_mpi = permute(mag2_mpi, ndims(mag2_mpi):-1:1);
        mag3_mpi = permute(mag3_mpi, ndims(mag3_mpi):-1:1);
        
        rho = squeeze(gas_mpi(1,:,:,:));
        vx = squeeze(gas_mpi(2,:,:,:));
        vy = squeeze(gas_mpi(3,:,:,:));
        vz = squeeze(gas_mpi(4,:,:,:));
        p = squeeze(gas_mpi(5,:,:,:));
        
        bi = mag1_mpi;
        bj = mag2_mpi;
        bk = mag3_mpi;
        
        I = 1:nx_total-1;
        J = 1:ny_total-1;
        K = 1:nz_total-1;
        
        b1 = zeros(size(x1ctr));
        b2 = zeros(size(x2ctr));
        b3 = zeros(size(x3ctr));
        B0_1 = zeros(size(x1ctr));
        B0_2 = zeros(size(x2ctr));
        B0_3 = zeros(size(x3ctr));
        
        b1(I,J,K) = bi(I,J,K) .* ((x1(I+1,J,K)-x1ctr(I,J,K)) ./ dx1(I,J,K)) + ...
                    bi(I+1,J,K) .* ((x1ctr(I,J,K)-x1(I,J,K)) ./ dx1(I,J,K));
        b2(I,J,K) = bj(I,J,K) .* ((x2(I,J+1,K)-x2ctr(I,J,K)) ./ dx2(I,J,K)) + ...
                    bj(I,J+1,K) .* ((x2ctr(I,J,K)-x2(I,J,K)) ./ dx2(I,J,K));
        b3(I,J,K) = bk(I,J,K) .* ((x3(I,J,K+1)-x3ctr(I,J,K)) ./ dx3(I,J,K)) + ...
                    bk(I,J,K+1) .* ((x3ctr(I,J,K)-x3(I,J,K)) ./ dx3(I,J,K));
        
        B0_1(I,J,K) = B01(I,J,K) .* ((x1(I+1,J,K)-x1ctr(I,J,K)) ./ dx1(I,J,K)) + ...
                      B01(I+1,J,K) .* ((x1ctr(I,J,K)-x1(I,J,K)) ./ dx1(I,J,K));
        B0_2(I,J,K) = B02(I,J,K) .* ((x2(I,J+1,K)-x2ctr(I,J,K)) ./ dx2(I,J,K)) + ...
                      B02(I,J+1,K) .* ((x2ctr(I,J,K)-x2(I,J,K)) ./ dx2(I,J,K));
        B0_3(I,J,K) = B03(I,J,K) .* ((x3(I,J,K+1)-x3ctr(I,J,K)) ./ dx3(I,J,K)) + ...
                      B03(I,J,K+1) .* ((x3ctr(I,J,K)-x3(I,J,K)) ./ dx3(I,J,K));
        
        rho_act = rho;
        v1_act = vx;
        v2_act = vy;
        v3_act = vz;
        p_act = p;
        
        b1_act = b1;
        b2_act = b2;
        b3_act = b3;
        
        B0_1_act = B0_1;
        B0_2_act = B0_2;
        B0_3_act = B0_3;
        
        bi_act = bi;
        bj_act = bj;
        bk_act = bk;
        
        h5_filename = sprintf('./data/data_%06d.h5', time_step);
        
        h5save(h5_filename, '/t', t);
        h5save(h5_filename, '/r', x1ctr_act);
        h5save(h5_filename, '/theta', x2ctr_act);
        h5save(h5_filename, '/phi', x3ctr_act);
        h5save(h5_filename, '/rho', rho_act);
        h5save(h5_filename, '/v1', v1_act);
        h5save(h5_filename, '/v2', v2_act);
        h5save(h5_filename, '/v3', v3_act);
        h5save(h5_filename, '/p', p_act);
        
        h5save(h5_filename, '/b1', b1_act);
        h5save(h5_filename, '/b2', b2_act);
        h5save(h5_filename, '/b3', b3_act);
        h5save(h5_filename, '/B01', B0_1_act);
        h5save(h5_filename, '/B02', B0_2_act);
        h5save(h5_filename, '/B03', B0_3_act);
        
        h5create(h5_filename, '/normalization_params', size(normalization_params));
        h5write(h5_filename, '/normalization_params', normalization_params);
        
        h5writeatt(h5_filename, '/normalization_params', 'description', ...
            'Normalization parameters: [x_Norm, u_Norm, Time_Norm, rho_Norm, p_Norm, B_Norm]');
        h5writeatt(h5_filename, '/normalization_params', 'units', ...
            '[m, m/s, s, kg/m^3, Pa, T]');
        
        fprintf('Successfully saved %s\n', h5_filename);
        
    catch ME
        fprintf('Error processing time step %d: %s\n', time_step, ME.message);
        continue; 
    end
end

fprintf('Processing completed! Processed %d time steps.\n', ww+1);
%%
function h5save(filename, dataset_path, data, compression_level)

    if nargin < 4
        compression_level = 6;  
    end
    
    data_size = size(data);
    
    chunk_size = calculate_chunk_size(data_size);
    
    h5create(filename, dataset_path, data_size, ...
        'Deflate', compression_level, ...
        'ChunkSize', chunk_size);
    
    h5write(filename, dataset_path, data);
end

function chunk_size = calculate_chunk_size(data_size)
    if numel(data_size) == 1
        chunk_size = min([1000, data_size]);
    elseif numel(data_size) == 2
        chunk_size = min([100, 100], data_size);
    elseif numel(data_size) == 3
        chunk_size = min([50, 50, 50], data_size);
    else
        chunk_size = min([50, 50, 50, ones(1, numel(data_size)-3)], data_size);
    end
end
%%
function reversed_array = reverseDimensions(array)
    ndims_array = ndims(array);
    order = ndims_array:-1:1;
    reversed_array = permute(array, order);
end
%%
function big_array = stitch_hdf5_array(num_i, num_j, num_k, crop_size, dataset_path,onface_i,onface_j,onface_k,start1,end1,seq)

    data_chunks = cell(num_i, num_j, num_k);

    for i = 0:num_i-1
        for j = 0:num_j-1
            for k = 0:num_k-1
                filename = sprintf('../build/mhd_%02d-%02d-%02d_%06d.h5', i, j, k, seq);

                data = h5read(filename, dataset_path);

                data_dims = ndims(data);

                
                switch data_dims
                    case 3
                        data  = data(:,:,:);
                        if k ~= 0
                            data = data(crop_size+1:end, :, :); 
                        end
                        if k ~= num_k-1
                            data = data(1:end-crop_size-1, :, :); 
                        end
                        if onface_k == 0 && k == num_k-1
                            data = data(1:end-1, :, :); 
                        end
                        
                        if j ~= 0
                            data = data(:, crop_size+1:end, :); 
                        end
                        if j ~= num_j-1
                            data = data(:, 1:end-crop_size-1, :); 
                        end
                        if onface_j == 0 && j == num_j-1
                            data = data(:,1:end-1, :); 
                        end
                        
                        if i ~= 0
                            data = data(:, :, crop_size+1:end); 
                        end
                        if i ~= num_i-1
                            data = data(:, :, 1:end-crop_size-1);
                        end
                        if onface_i == 0 && i == num_i-1
                            data = data(:,:,1:end-1); 
                        end
                        
                    case 4
                        data  = data(:,:,:,start1:end1);
                        if k ~= 0
                            data = data(crop_size+1:end, :, :, :); 
                        end
                        if k ~= num_k-1
                            data = data(1:end-crop_size-1, :, :, :); 
                        end
                        if onface_k == 0 && k == num_k-1
                            data = data(1:end-1, :, :, :); 
                        end
                        
                        if j ~= 0
                            data = data(:, crop_size+1:end, :, :); 
                        end
                        if j ~= num_j-1
                            data = data(:, 1:end-crop_size-1, :, :); 
                        end
                        if onface_j == 0 && j == num_j-1
                            data = data(:,1:end-1, :, :); 
                        end
                        
                        if i ~= 0
                            data = data(:, :, crop_size+1:end, :); 
                        end
                        if i ~= num_i-1
                            data = data(:, :, 1:end-crop_size-1, :);
                        end
                        if onface_i == 0 && i == num_i-1
                            data = data(:,:,1:end-1, :); 
                        end
                        
                    case 5
                        data  = data(:,:,:,:,start1:end1);
                        if k ~= 0
                            data = data(crop_size+1:end, :, :, :,:); 
                        end
                        if k ~= num_k-1
                            data = data(1:end-crop_size-1, :, :, :,:); 
                        end
                        if onface_k == 0 && k == num_k-1
                            data = data(1:end-1, :, :, :,:); 
                        end
                        if j ~= 0
                            data = data(:, crop_size+1:end, :, :,:); 
                        end
                        if j ~= num_j-1
                            data = data(:,1:end-crop_size-1, :, :,:); 
                        end
                        if onface_j == 0 && j == num_j-1
                            data = data(:,1:end-1, :, :,:); 
                        end
                        if i ~= 0
                            data = data(:, :, crop_size+1:end, :,:); 
                        end
                        if i ~= num_i-1
                            data = data(:, :, 1:end-crop_size-1, :,:); 
                        end
                        if onface_i == 0 && i == num_i-1
                            data = data(:,:,1:end-1, :,:); 
                        end
                        
                    otherwise
                        error('Unsupported data dimensionality: %d. Expected 4D or 5D or 3D data.', data_dims);
                end
                data_chunks{i+1, j+1, k+1} = data;
            end
        end
    end

    concat_i = cell(num_j, num_k); 

    for j = 1:num_j
        for k = 1:num_k
            concat_i{j, k} = cat(3, data_chunks{1:num_i, j, k}); 
        end
    end

    concat_j = cell(1, num_k);
    for k = 1:num_k
        concat_j{k} = cat(2, concat_i{1:num_j, k});
    end

    big_array = cat(1, concat_j{:}); 
end
%%
function max_seq = find_max_seq()
    files = dir('../build/mhd_*-*-*_*.h5');
    
    max_seq = -Inf;
    
    for k = 1:length(files)
        filename = files(k).name;
        
        tokens = regexp(filename, '^mhd_\d{2}-\d{2}-\d{2}_(\d{6})\.h5$', 'tokens');
        
        if ~isempty(tokens)
            seq = str2double(tokens{1}{1});
            
            if seq > max_seq
                max_seq = seq;
            end
        end
    end

    if max_seq == -Inf
        error('There is no files！');
    end
end
