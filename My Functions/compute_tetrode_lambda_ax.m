function [tetrode_lamb_ax, X_bin, Y_bin, Z_bin] = compute_tetrode_lambda_ax(bw, encoding_vectors, mu, pixels_per_m1, n_bins_y, n_bins_x, xy1, n_bins_z, pi_x, bin_size1, n_tetrodes, n_tetrode_spikes)
tetrode_lamb_ax = zeros(n_bins_y, n_bins_x, n_bins_z, size(n_tetrode_spikes, 1));
tic
for tetrode = 1:n_tetrodes
    disp(tetrode)
    if tetrode == 1
        [lamb_ax, X_bin, Y_bin, Z_bin] = lambda_ax_precomputed(bw, encoding_vectors(1:n_tetrode_spikes(1), :), mu, pixels_per_m1, n_bins_y, n_bins_x, xy1, n_bins_z, pi_x, bin_size1);
    else
        [lamb_ax, X_bin, Y_bin, Z_bin] = lambda_ax_precomputed(bw, encoding_vectors(sum(n_tetrode_spikes(1:tetrode-1)):sum(n_tetrode_spikes(1:tetrode)), :), mu, pixels_per_m1, n_bins_y, n_bins_x, xy1, n_bins_z, pi_x, bin_size1);
    end
    tetrode_lamb_ax(:, :, :, tetrode) = lamb_ax;
end
toc
end