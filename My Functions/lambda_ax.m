function lamb_ax = lambda_ax(X_bin_grid, Y_bin_grid, bw, enc_spike_vectors, mu, pixels_per_m1, max_time, n_bins_y, n_bins_x, xy1, dec_spike_amps, pi_x)

X_bin_vec = X_bin_grid(:); % It flattens the array column wise
Y_bin_vec = Y_bin_grid(:);
XY_bins_vec = [X_bin_vec, Y_bin_vec]; % Two column vectors of all coordinates
A_mult_vec = repmat(dec_spike_amps, size(XY_bins_vec, 1), 1); % so its repeated for each bin
AXY_bins_vec = [A_mult_vec XY_bins_vec]; % Locations to be evaluated
pi_x = repmat(pi_x, size(AXY_bins_vec, 1)/size(pi_x(:), 1), 1);
p_x  = mvksdensity(enc_spike_vectors, AXY_bins_vec, 'Bandwidth', bw); % KS density of firing rates in specific bin
%pi_x = mvksdensity(xy1(1:max_time, :)/pixels_per_m1*100, XY_bins_vec, 'Bandwidth', bw(); % KS density of position map
disp(size(p_x));
disp(size(pi_x));
lamb_ax = mu * p_x ./ pi_x; % Compute the firing rate map that only depends on position
lamb_ax = reshape(lamb_ax, n_bins_y, n_bins_x); % Reshape so its a meshgrid for plotting
