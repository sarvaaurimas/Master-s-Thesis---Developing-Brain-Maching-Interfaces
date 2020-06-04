function [log_post, decoded_position, n_enc_spikes, log_sum_lamb_ax_arr] = compute_unsorted_posterior_position_precomp(decoding_data, dt, n_bins_x, n_bins_y, n_tetrodes, n_electrodes, pixels_per_m1, lamb_x, bin_size1, bin_z_coord, tetrode_lamb_ax, X_bin, Y_bin, interspike_interval)
samp_rate = 30000;
n_time_bins = floor(size(decoding_data, 2)/samp_rate/dt); % round down to compute number of time bins
log_post = zeros(n_bins_y, n_bins_x, n_time_bins); % Empty array to hold the posterior
decoded_position = zeros(n_time_bins, 2);

n_enc_spikes = zeros(n_time_bins, n_tetrodes);
log_sum_lamb_ax_arr = zeros(n_bins_y, n_bins_x, n_time_bins);


for t_bin = 1:n_time_bins
    disp(t_bin);
    %disp('% done');
    %disp(t_bin/n_time_bins*100);
    % Compute where each bin starts and ends in the sampling space
    t_bin_start = 1 + (t_bin - 1) * samp_rate * dt;
    t_bin_end = t_bin * samp_rate * dt;
    
    log_bin_probability = zeros(n_bins_y, n_bins_x);
    
    % Placeholder to count spikes in previous tetrodes
    %prev_spikes = 0;
    
    for tetrode = 1:n_tetrodes
        disp('Tetrode');
        disp(tetrode);
        % only electrodes belocking to that tetrode in that specific time
        % bin
        tetrode_data = decoding_data(1 + (tetrode - 1) * 4 : tetrode * 4,  t_bin_start:t_bin_end);
        tetrode_data = tetrode_data(1:n_electrodes, :);
        % Compute tetrode spike times 
        tetrode_spike_time = spike_times_v2(tetrode_data, interspike_interval); % obtain the spike times
        %disp(size(tetrode_spike_time));
        % tetrode amplitudes are 4xn matrix where n is the spike times
        tetrode_amplitudes = tetrode_data(:, tetrode_spike_time)';
        n_spikes = size(tetrode_spike_time, 2);
        
        
        % Choose only the tetrode encoding vectors
        %curr_tetrode_spikes = n_tetrode_spikes(tetrode);
        %tetrode_encoding_vectors = encoding_vectors((1 + prev_spikes):(prev_spikes + curr_tetrode_spikes), :);
        %prev_spikes = prev_spikes + curr_tetrode_spikes;
        
        %% The following part computes sum of lamb ax for each spike
        sum_log_lamb_ax = zeros(n_bins_y, n_bins_x);
        
        
        
        for spike_ind = 1:n_spikes
           %disp(spike_ind);
           % The electode amplitudes for one spike 1xn_features vector
           current_amplitude = roundtowardvec(tetrode_amplitudes(spike_ind, 1:n_electrodes), bin_z_coord);
           disp([tetrode_amplitudes(spike_ind, 1:n_electrodes), current_amplitude]);
           % Compute lambda ax map for each spike
           lamb_ax = tetrode_lamb_ax(:, :, bin_z_coord == current_amplitude, tetrode);
           sum_log_lamb_ax = sum_log_lamb_ax + log(lamb_ax);
           %surf(X_bin(:, :, 1), Y_bin(:, :, 1), log(lamb_ax));
           %pause;
        end
        % Update the probability
        log_prob_curr_tetrode = n_spikes * log(dt) + sum_log_lamb_ax - dt * lamb_x;
        % Add current tetrode probability as in eq 4 in Cloosterman
        log_bin_probability = log_bin_probability + log_prob_curr_tetrode; 
        % Update n spikes array
        n_enc_spikes(t_bin, tetrode) = n_spikes;
        
        %surf(X_bin(:, :, 1), Y_bin(:, :, 1), sum_log_lamb_ax);
        %pause;
        surf(X_bin(:, :, 1), Y_bin(:, :, 1), log_bin_probability);
        pause;
        
    end
    
    % Update the log sum ax array
    log_sum_lamb_ax_arr(:, :, t_bin) = sum_log_lamb_ax;
    %surf(X_bin(:, :, 1), Y_bin(:, :, 1), sum_log_lamb_ax);
    %disp(sum_log_lamb_ax);
    %pause;
    
    log_post(:, :, t_bin) = log_bin_probability;
    
    max_bin_post = max(max(log_bin_probability)); % Calculate maximum posteriors in each step
    [row, col] = find(log_bin_probability == max_bin_post); % Get indexes of max posteriors in each step, sometimes throws errors if using too little cells
   
    y_bin = (row-0.5)*bin_size1/pixels_per_m1*100;
    x_bin = (col-0.5)*bin_size1/pixels_per_m1*100; % Predicted centers of bins in pixels
    
    decoded_position(t_bin, 1) = x_bin;
    decoded_position(t_bin, 2) = y_bin;
end

end