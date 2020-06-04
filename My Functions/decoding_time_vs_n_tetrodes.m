function times = decoding_time_vs_n_tetrodes(n_electrodes, interspike_interval, dec_raw_data, pixels_per_m1, X_bin, Y_bin, n_bins_y, n_bins_x, bin_size1, dt, tetrode_lamb_x, tetrode_z_coord, tetrode_lamb_ax, sigmas, n_step_vel)
tetrodes_n = [2:16];
times = zeros(1, size(tetrodes_n, 2));
for i = 1:size(tetrodes_n, 2)
    tic
    n_tetrodes = tetrodes_n(i);
    
    [log_post, decoded_position, n_dec_spikes, log_sum_lamb_ax_arr, dec_spikes] = compute_unsorted_posterior_position_precomp_v3(...
                                                           dec_raw_data, dt, n_bins_x, n_bins_y, n_tetrodes, ...
                                                       n_electrodes, pixels_per_m1, tetrode_lamb_x, bin_size1, ...
                                                       tetrode_z_coord, tetrode_lamb_ax, ...
                                                        interspike_interval, sigmas, X_bin, Y_bin, n_step_vel);
    times(i) = toc  
end

end