function plot_ax_slice_z_var(tetrode_lamb_ax, X_bin, Y_bin, Z_bin)

for tetrode = 1:size(tetrode_lamb_ax, 4)
    disp('tetrode');
    disp(tetrode);
    for slice = 1:size(tetrode_lamb_ax, 3)
        disp(slice);
        disp(size(tetrode_lamb_ax));
        lamb_ax = reshape(tetrode_lamb_ax(:, :, :, tetrode), size(tetrode_lamb_ax, 1), size(tetrode_lamb_ax, 2), size(tetrode_lamb_ax, 3));
        X = reshape(X_bin(:, :, slice), size(X_bin, 1),size(X_bin, 2));
        Y = reshape(Y_bin(:, :, slice), size(Y_bin, 1),size(Y_bin, 2));
        lambd = reshape(lamb_ax(:, :, slice), size(lamb_ax, 1),size(lamb_ax, 2));
        surf(X, Y, lambd);
        xlabel('X, cm');
        ylabel('Y, cm');
        zlabel('Firing rate, Hz');
        title(Z_bin(1, 1, slice));
        pause;
    end
end

end