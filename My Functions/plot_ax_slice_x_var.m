function plot_ax_slice_x_var(tetrode_lamb_ax, X_bin, Y_bin, Z_bin, n_slice)

for tetrode = 1:size(tetrode_lamb_ax, 4)
    lamb_ax = reshape(tetrode_lamb_ax(:, :, :, tetrode), size(tetrode_lamb_ax, 1), size(tetrode_lamb_ax, 2), size(tetrode_lamb_ax, 3));
    X = reshape(X_bin(n_slice, :, :), size(X_bin, 2),size(X_bin, 3));
    Y = reshape(Z_bin(n_slice, :, :), size(Z_bin, 2),size(Z_bin, 3));
    lambd = reshape(lamb_ax(n_slice, :, :), size(lamb_ax, 2),size(lamb_ax, 3));
    size(X);
    size(Y);
    size(lambd);
    surf(X, Y, lambd);
    %surf(reshape(X_bin(n_slice, :, :), size(X_bin, 2),size(X_bin, 3)) , reshape(Z_bin(n_slice, :, :), size(Z_bin, 2),size(Z_bin, 3)), reshape(lamb_ax(n_slice, :, :), size(lamb_ax, 2),size(lamb_ax, 3)));
    xlabel('X, cm');
    ylabel('Amplitude, \mu V');
    zlabel('Firing rate, Hz');
    pause;
end
end