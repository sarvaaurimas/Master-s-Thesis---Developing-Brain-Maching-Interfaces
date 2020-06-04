function plot_tetrodes(t)

% Plots the spike reading for t timestamps and k channels
samp_freq = 30000;
t = t*samp_freq;
fid = fopen('spikes.bin','r');
dat_t = (fread(fid,[1, 64*t], '*int16'));
data=reshape(dat_t,64,length(dat_t)/64); 
data_offset=double(data)+double(repmat([(0:63)*200]',1,size(data,2)));
data_offset = double(data_offset);

time = [1:size(data_offset, 2)]/samp_freq * 1000; % in ms


for tetrode = 1:16
    disp(tetrode);
    plot(time, data_offset((tetrode-1)*4 + 1:(tetrode)*4, :)');
    xlabel('Time, ms','FontSize',14);
    ylabel('Voltage, \mu V','FontSize',14);
    pause();
end

