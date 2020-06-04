function tetrode_lamb_x = tetrode_lambda_x(X_bin_grid, Y_bin_grid, bw, n_bins_y, n_bins_x, encoding_vectors, n_tetrode_spikes, pi_x, n_tetrodes, tetrode_mu)

X_bin_vec = X_bin_grid(:); % It flattens the array column wise
Y_bin_vec = Y_bin_grid(:);
XY_bins_vec = [X_bin_vec, Y_bin_vec]; % Two column vectors of all coordinates

tetrode_lamb_x = zeros(n_bins_y, n_bins_x, n_tetrodes);
for tetrode = 1:n_tetrodes
    disp('Tetrode');
    disp(tetrode);
    % only electrodes belocking to that tetrode in that specific time
    % bin
    tic
    if tetrode == 1
        p_x  = mvksdensity(encoding_vectors(1:n_tetrode_spikes(1), 2:3), XY_bins_vec, 'Bandwidth', bw); % KS density of firing rates in specific bin
    else
        p_x  = mvksdensity(encoding_vectors(sum(n_tetrode_spikes(1:tetrode-1)):sum(n_tetrode_spikes(1:tetrode)), 2:3), XY_bins_vec, 'Bandwidth', bw); % KS density of firing rates in specific bin
    end
    
    lamb_x = tetrode_mu(tetrode) * p_x ./ pi_x; % Compute the firing rate map that only depends on position

    lamb_x = reshape(lamb_x, n_bins_y, n_bins_x); % Reshape so its a meshgrid for plotting
    tetrode_lamb_x(:, :, tetrode) = lamb_x;
    toc

end
