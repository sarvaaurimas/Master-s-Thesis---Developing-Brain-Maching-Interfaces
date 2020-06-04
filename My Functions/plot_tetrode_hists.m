function plot_channel_hists(t)

% Plots the spike reading for t timestamps and k channels

fid = fopen('spikes.bin','r');
dat_t = (fread(fid,[1, 64*t], '*int16'));
data=reshape(dat_t,64,length(dat_t)/64); 
data_offset=double(data)+double(repmat([(0:63)*200]',1,size(data,2)));
data_offset = double(data_offset);
data = double(data);
%for tetrode = 1:16
    %disp(tetrode);
    %plot(data_offset((tetrode-1)*4 + 1:(tetrode)*4, :)');
    %pause();
%end

T = [1:t]/30000; % the time array for plotting

for channel = 1:64
    
    
    disp(channel);
    subplot(2,1,1);
    histfit(data(channel, :), 20);
    sigma = std(data(channel, :));
    title(['Standard deviation = ', num2str(sigma), '\muV']);
    xlabel('Channel potential, \muV','FontSize',12);
    ylabel('Number of occurances','FontSize',12);
    legend('Channel voltage distribution','Gaussian fit');
    
    
    subplot(2,1,2); 
    sigma_arr = ones(1, size(data_offset, 2))*4*sigma;
    plot(T, data(channel, :)');
     
    disp(size(sigma_arr));
    hold on;
    plot(T, sigma_arr);
    hold on;
    plot(T, -sigma_arr);
    hold off;
    legend('Channel voltage','+4 standard deviations', '-4 standard deviations');
    xlabel('Time, s','FontSize',12);
    ylabel('Channel potential, \muV','FontSize',12);
    pause();
end