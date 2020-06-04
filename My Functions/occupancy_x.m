function pi_x = occupancy_x(X_bin_grid, Y_bin_grid, bw, pixels_per_m1, training_time, n_bins_y, n_bins_x, xy1)

X_bin_vec = X_bin_grid(:); % It flattens the array column wise
Y_bin_vec = Y_bin_grid(:);
XY_bins_vec = [X_bin_vec, Y_bin_vec]; % Two column vectors of all coordinates

pi_x = mvksdensity(xy1(1:training_time*50, :)/pixels_per_m1*100, XY_bins_vec, 'Bandwidth', bw); % KS density of position map

end