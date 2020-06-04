function spikes = load_spikes(t)

t = 30000*t; % as its sampled at 30 kHz and t is in seconds
fid = fopen('spikes.bin','r');
dat_t = (fread(fid,[1, 64*t], '*int16'));
data=reshape(dat_t,64,length(dat_t)/64); 
spikes = double(data);