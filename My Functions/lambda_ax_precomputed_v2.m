function [lamb_ax, X_bin, Y_bin, Z_bin, bin_z_coord] = lambda_ax_precomputed_v2(bw, enc_spike_vectors, mu, pixels_per_m1, n_bins_y, n_bins_x, n_bins_z, pi_x, bin_size1)

bin_x_coord = ([1 : n_bins_x] -0.5)*bin_size1/(pixels_per_m1/100); % Coodinates of each bin
bin_y_coord = ([1 : n_bins_y] -0.5)*bin_size1/(pixels_per_m1/100);

z_min = min(enc_spike_vectors(:, 1));
z_max = max(enc_spike_vectors(:, 1));
disp(z_min);
disp(z_max);
bin_z_coord = linspace(z_min, z_max, n_bins_z); % From min to max values
[X_bin, Y_bin, Z_bin] = meshgrid(bin_x_coord, bin_y_coord, bin_z_coord); % A meshgrid for each bin center




X_bin_vec = X_bin(:); % It flattens the array column wise
Y_bin_vec = Y_bin(:);
Z_bin_vec = Z_bin(:);
AXY_bins_vec = [Z_bin_vec, X_bin_vec, Y_bin_vec]; % Two column vectors of all coordinates


disp(size(enc_spike_vectors));
disp(size(AXY_bins_vec));
p_x  = mvksdensity(enc_spike_vectors, AXY_bins_vec, 'Bandwidth', bw); % KS density of firing rates in specific bin

p_x = reshape(p_x, n_bins_y, n_bins_x, n_bins_z); % Reshape so its a meshgrid for plotting
pi_x = repmat(reshape(pi_x, n_bins_y, n_bins_x), 1, 1, n_bins_z);

lamb_ax = mu * p_x ./ pi_x; % Compute the firing rate map that only depends on position
