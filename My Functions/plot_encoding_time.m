function times = plot_encoding_time(n_tetrodes, n_electrodes, interspike_interval, enc_data, xy1, pixels_per_m1, B_dec, B_x, X_bin, Y_bin, n_bins_y, n_bins_x, n_bins_z, bin_size1)
encoding_minutes = [1:9]
times = zeros(1, size(encoding_minutes, 2));
for i = 1:size(encoding_minutes, 2);
    tic
    curr_encoding_time = encoding_minutes(i)*60;
    curr_encoding_data = enc_data(:, 1:curr_encoding_time*30000);
    [encoding_vectors, n_tetrode_spikes, sigmas] = compute_encoding_vectors_v2(curr_encoding_data, n_tetrodes, n_electrodes, curr_encoding_time, xy1, pixels_per_m1, interspike_interval);

    training_pi_time = curr_encoding_time;


    pi_x = occupancy_x(X_bin, Y_bin, B_x, pixels_per_m1, training_pi_time, n_bins_y, n_bins_x, xy1);

    tetrode_mu = n_tetrode_spikes/curr_encoding_time; % Array of inidividual tetrode mu's
    tetrode_lamb_x = tetrode_lambda_x(X_bin, Y_bin, B_x, n_bins_y, n_bins_x, encoding_vectors, n_tetrode_spikes, pi_x, n_tetrodes, tetrode_mu);
    [tetrode_lamb_ax, X_bin_3d, Y_bin_3d, Z_bin_3d, tetrode_z_coord] = compute_tetrode_lambda_ax_v2(B_dec, encoding_vectors, tetrode_mu, pixels_per_m1, n_bins_y, n_bins_x, n_bins_z, pi_x, bin_size1, n_tetrodes, n_tetrode_spikes);
    times(i) = toc  
end

end