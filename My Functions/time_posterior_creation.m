function [cells, elapsed_time] = time_posterior_creation(n_test_bins, obsSpk, ratemaps, posmap, tBin, bin_size1, pixels_per_m1, n_cells)
n = 10;
cells = linspace(1, 255, n);
elapsed_time = zeros(1,n);

for i = 1:n
   disp(i);
   cell_n = cells(i);
   rms = ratemaps(:, :, 1:cell_n); 
   tic 
   [x, y] = decode_posPred(n_test_bins, obsSpk, rms, posmap, tBin, bin_size1, pixels_per_m1, cell_n);
   elapsed_time(i) = toc;
end
end