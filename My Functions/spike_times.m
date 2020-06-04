function [unique_spike_times, difference] = spike_times(spike_data, sim_spike_distance)
spike_data = -spike_data; % negative as we want to find peaks not bottoms
sptimes = [];
spchannels = [];

for channel = 1:size(spike_data, 1)
   % for each channel compute peak locations
   sigma = std(spike_data(channel, :));
   [~,chan_sp_times] = findpeaks(spike_data(channel, :),'MinPeakHeight',4*sigma); % Use 3.5 standard deviations as the threshold
   sptimes = cat(2, sptimes, chan_sp_times);
   channel_times = ones(1, size(chan_sp_times, 2))*channel;
   spchannels = cat(2, spchannels, channel_times);
end

[sptimes, indx] = sort(sptimes); % return all unique times of thresholded spikes
spchannels = spchannels(indx); % sort the channels array too
[sptimes, ia, ic] = unique(sptimes); % get unique times
spchannels = spchannels(ia); % channels of each spike

% THIS PART CALCULATES UNIQUE SPIKE TIMES
unique_spike_times = [];
i = 1;
difference = [];
while_end = 0;

while i<size(sptimes, 2) % Run it until the end of the spike times vector
    diff = sptimes(i+1) - sptimes(i); % Compute the difference 
    difference = cat(2, difference, [diff]);
    
    if while_end ==1
        while_end = 0; % dont add it to the unique spike times
        sim_spike_times = [];
    else
       sim_spike_times = [sptimes(i)]; % store the timestamps of peaks belonging to the same spike here
    end
    
    
    while diff < sim_spike_distance && i < (size(sptimes, 2)-1)
        sim_spike_times = cat(2, sim_spike_times, [sptimes(i+1)]); % store the new peak timestamp here
        i = i+1;
        diff = sptimes(i+1) - sptimes(i); % Compute the difference
        while_end = 1;
    end
    
    
    if while_end == 1
        i = i-1;
        sim_spike_times = cat(2, sim_spike_times, [sptimes(i+1)]); %add the last one that otherwise would not be included in the loop
    end
    
    
    unique_spike_times = cat(2, unique_spike_times, [round(mean(sim_spike_times))]); % Add the mean of the peaks as the new time
    i = i+1;
end
unique_spike_times = unique_spike_times(~isnan(unique_spike_times));



