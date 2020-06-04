function [x, y, row, col, post, max_post] = decode_posPred(n_test_bins, obsSpk, ratemaps, posmaps, tBin, bin_size1, pixels_per_m1, n_cells, speedScaleFact)
if nargin < 9
    speedScaleFact = 1;
end

post = decode_calcBayesPost(obsSpk(1:n_cells, 1:n_test_bins), ratemaps(:, :, 1:n_cells), posmaps, tBin, speedScaleFact); % Compute posteriors for each step 
max_post = reshape(max(max(post)), [n_test_bins, 1]); % Calculate maximum posteriors in each step
row = zeros(length(max_post), 1);
col = zeros(length(max_post), 1);

for i = 1:length(max_post)
    [row(i), col(i)] = find(post(:, :, i) == max_post(i)); % Get indexes of max posteriors in each step, sometimes throws errors if using too little cells
end
y = (row-0.5)*bin_size1/pixels_per_m1*100;
x = (col-0.5)*bin_size1/pixels_per_m1*100; % Predicted centers of bins in pixels
