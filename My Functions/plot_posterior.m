function F = plot_posterior(post)
% frame = 141;
% curr_frame = post(:, :, frame);
% curr_frame(curr_frame < 0.001) = NaN;
% curr_frame(isnan(curr_frame))=0;
x = [1:132];
y = [1:113];
[X, Y] = meshgrid(x, y);
% surf(X, Y, curr_frame, 'Facecolor', 'interp');
% colorbar;
% zlim([0.01, max(max(curr_frame))]);


%Z = peaks;
%surf(Z)
%axis tight manual
%ax = gca;
%ax.NextPlot = 'replaceChildren';


loops = 40;
F(loops) = struct('cdata',[],'colormap',[]);
for frame = 1:loops
    curr_frame = post(:, :, frame);
    curr_frame(curr_frame < 0.001) = NaN;
    curr_frame(isnan(curr_frame))=0;
    surf(X,Y, curr_frame, 'Facecolor', 'interp')
    zlim([0.01, 0.8]);
    xlabel('x, cm','FontSize',14);
    ylabel('y, cm','FontSize',14);
    zlabel('Posterior probability','FontSize',14);
    %xticklabels([0, 20, 40, 60, 80, 100, 120]*2.5);
    %yticklabels([1:11]*10*2.5)
    
    drawnow
    F(frame) = getframe;
    pause;

end
