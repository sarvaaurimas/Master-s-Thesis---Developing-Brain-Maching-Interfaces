function [lamb_x, pi_x] = lambda_x(X_bin_grid, Y_bin_grid, bw, spike_coords, mu, pixels_per_m1, max_time, n_bins_y, n_bins_x, xy1)

X_bin_vec = X_bin_grid(:); % It flattens the array column wise
Y_bin_vec = Y_bin_grid(:);
XY_bins_vec = [X_bin_vec, Y_bin_vec]; % Two column vectors of all coordinates

p_x  = mvksdensity(spike_coords, XY_bins_vec, 'Bandwidth', bw); % KS density of firing rates in specific bin
pi_x = mvksdensity(xy1(1:max_time, :)/pixels_per_m1*100, XY_bins_vec, 'Bandwidth', bw); % KS density of position map

lamb_x = mu * p_x ./ pi_x; % Compute the firing rate map that only depends on position
lamb_x = reshape(lamb_x, n_bins_y, n_bins_x); % Reshape so its a meshgrid for plotting


