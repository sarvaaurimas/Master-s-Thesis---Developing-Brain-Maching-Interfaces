function data = plot_channels(t, channel_start, channel_end, sptime)

% Plots the spike reading for t timestamps and k channels

fid = fopen('spikes.bin','r');
dat_t = (fread(fid,[1, 64*t], '*int16'));
data=reshape(dat_t,64,length(dat_t)/64); 
data_offset=double(data)+double(repmat([(0:63)*200]',1,size(data,2)));
%n_spikes = size(sptime < t, 1); % how many spikes occur in t period
%det_spikes = zeros(64, t); 
%det_spikes(:, sptime(sptime<t)) = 100; %set the detected spikes matrix to be equal to 100 mv on spike
%det_spikes =det_spikes+double(repmat([(0:63)*200]',1,size(det_spikes,2)));
data_offset = double(data_offset);
data = double(data);
%det_spikes = double(det_spikes);
plot(data_offset(channel_start:channel_end, :)');
%hold on;
%plot(det_spikes(channel_start:channel_end, :)');