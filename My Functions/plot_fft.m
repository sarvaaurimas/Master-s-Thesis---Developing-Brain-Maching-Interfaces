function plot_fft(spikes)
Fs = 30000; % 30 Khz sampling frequency
L = size(spikes, 2);
amps = fft(spikes);
amps = abs(amps/L);
P1 = amps(1:L/2+1);
f = Fs*(0:(L/2))/L;
plot(f, P1);