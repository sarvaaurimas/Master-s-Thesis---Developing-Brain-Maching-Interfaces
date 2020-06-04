function [encoding_vectors, tetrode_spike_n, sigmas, enc_spike_n] = compute_encoding_vectors_v2(spike_data, n_tetrodes, n_electrodes, enc_time, position, pixels_per_m1, interspike_interval, dt, x_real, y_real)
samp_rate = 30000;
enc_time = enc_time * samp_rate;
tetrode_spike_n = zeros(n_tetrodes, 1);
encoding_vectors = zeros(0, n_electrodes + 2); % 2 position coordinates
sigmas = zeros(n_tetrodes, 1);
for tetrode = 1:n_tetrodes
    disp('Tetrode');
    disp(tetrode);
    
    tetrode_data = spike_data((1 + (tetrode - 1) * 4) : (tetrode * 4),  1:enc_time);
    tetrode_data = tetrode_data(1:n_electrodes, :);
    
    %(size(tetrode_data));
    [tetrode_spike_time, sigma] = spike_times_v2(tetrode_data, interspike_interval);  
    sigmas(tetrode) = sigma;
    
    % Compute the tetrode amplitudes, n x n_electrode vector
    tetrode_amplitudes = tetrode_data(:, tetrode_spike_time)';
    
    % Compute how many spikes per tetrode
    n_spikes = size(tetrode_spike_time, 2);
    
    % Assign the tetrode spikes to count array
    tetrode_spike_n(tetrode) = n_spikes;
    
    % Get the position of each spike, n x 2 vector
    spike_positions = position(ceil(tetrode_spike_time/samp_rate*50), :);
    
    tetrode_encoding_vectors = [tetrode_amplitudes spike_positions/pixels_per_m1*100];
    encoding_vectors = vertcat(encoding_vectors, tetrode_encoding_vectors);
    
       
end
end