function [log_post, decoded_position] = compute_unsorted_posterior_precomp(spike_data, dt, lamb_x, encoding_vectors, n_bins_x, n_bins_y, n_tetrodes, spike_distance, n_of_electrodes, X_bin, Y_bin, bw, mu, pixels_per_m1, train_time, bin_z_coord)
samp_rate = 30000;
n_time_bins = floor(size(spike_data, 2)/samp_rate/dt); % round down to compute number of time bins
log_post = zeros(n_bins_y, n_bins_x, n_time_bins); % Empty array to hold the posterior
decoded_position = zeros(n_time_bins, 2);



for t_bin = 1:n_time_bines
    
    % Compute where each bin starts and ends in the sampling space
    t_bin_start = 1 + (t_bin - 1) * samp_rate * dt;
    t_bin_end = t_bin * samp_rate * dt;
    
    log_bin_probability = zeros(n_bins_y, n_bins_x);
    for tetrode = 1:n_tetrodes
        disp('Tetrode');
        disp(tetrode);
        % only electrodes belocking to that tetrode in that specific time
        % bin
        tetrode_data = spike_data(1 + (tetrode - 1) * 4 : tetrode * 4,  t_bin_start:t_bin_end);
        disp(size(tetrode_data));
        % Compute tetrode spike times 
        [tetrode_spike_time, difference] = spike_times(tetrode_data, spike_distance); % obtain the spike times
        % tetrode amplitudes are 4xn matrix where n is the spike times
       
        
        % Choose only the tetrode encoding vectors
        tetrode_encoding_vectors = encoding_vectors(:, :, tetrode);
        %% The following part computes sum of lamb ax for each spike
        sum_lamb_ax = zeros(n_bins_y, n_bins_x);
        for spike_ind = 1:n_spikes
           % The electode amplitudes for one spike 1xn_features vector
           current_amplitudes = tetrode_amplitudes(1:n_of_electrodes , spike_ind)';
           % Compute lambda ax map for each spike
           lamb_ax = lambda_ax(X_bin, Y_bin, bw, encoding_vectors, mu, pixels_per_m1, train_time, n_bins_y, n_bins_x, xy1, current_amplitudes);
           sum_lamb_ax = sum_lamb_ax + lamb_ax
        end
        log_prob_curr_tetrode = n_spikes * log(dt) + sum_lamb_ax - dt * lamb_x;
        % Add current tetrode probability as in eq 4 in Cloosterman
        log_bin_probability = log_bin_probability + log_prob_curr_tetrode;            
    end
    
    log_post(:, :, t_bin) = log_bin_probability
    
    max_bin_post = max(max(log_bin_probability)); % Calculate maximum posteriors in each step
    [row, col] = find(log_bin_probability == max_bin_post); % Get indexes of max posteriors in each step, sometimes throws errors if using too little cells
    y_bin = (row-0.5)*bin_size1/pixels_per_m1*100;
    x_bin = (col-0.5)*bin_size1/pixels_per_m1*100; % Predicted centers of bins in pixels
    
    decoded_position(t_bin, 1) = x_bin;
    decoded_position(t_bin, 2) = y_bin;
end

end