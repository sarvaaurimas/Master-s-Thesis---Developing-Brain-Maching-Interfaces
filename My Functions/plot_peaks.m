function plot_peaks(spike_data, threshold)

data_offset=spike_data+double(repmat([(0:size(spike_data, 1))*200]',1,size(spike_data,2)));


[~,locs_toppeak] = findpeaks(spike_data,'MinPeakHeight',threshold);
spikes_inverted = -spike_data;
[~,locs_botpeak] = findpeaks(spikes_inverted,'MinPeakHeight',threshold);
disp(size(data_offset));
figure
hold on 
plot(data_offset')
plot(locs_toppeak,spike_data(:, locs_toppeak)','rv','MarkerFaceColor','r')
plot(locs_botpeak,spike_data(:, locs_botpeak)','rs','MarkerFaceColor','b')
%axis([0 1850 -1.1 1.1])
grid on
legend('Original signal','Top Peaks','Bottom peaks')
xlabel('Samples')
ylabel('Voltage(mV)')
title('R-wave and S-wave in Noisy ECG Signal')