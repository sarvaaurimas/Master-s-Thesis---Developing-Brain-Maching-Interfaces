function [unique_spike_times, sigma] = spike_times_v2(spike_data, spike_latency_interval, sigma)
spike_data = -spike_data; % negative as we want to find peaks not bottoms
unique_spike_times = [];

if nargin < 3
    sigma = std(spike_data);
end
%disp(4*sigma);
[~, sptimes] = findpeaks(spike_data,'MinPeakHeight',4*sigma); % Use 3.5 standard deviations as the threshold
i = 1;
interspike_interval = ceil(spike_latency_interval / 1000 * 30000); % latency interval is in ms
prev_false = 0;
while i<size(sptimes, 2)
    
    diff = sptimes(i+1) - sptimes(i);
    if diff <= interspike_interval
        i = i+1;
        prev_false = 1;
    else
        if not(prev_false)
           unique_spike_times = cat(2, unique_spike_times, [sptimes(i)]);
        end
        i = i+1;
        prev_false = 0;
    end
     
end
end
