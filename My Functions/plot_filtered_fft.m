function filtered_spikes = filtered_fft(spikes, bands)
Fs = 30000; % 30 Khz sampling frequency
%L = size(spikes, 2);
%amps = fft(spikes);
%amps = abs(amps/L);
%P1 = amps(1:L/2+1);
filtered_spikes = highpass(spikes, bands, Fs);
%f = Fs*(0:(L/2))/L;
%plot(f, P1);
