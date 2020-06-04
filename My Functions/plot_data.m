fid = fopen('spikes.bin','r');
dat_t = fread(fid,'*int16');
data=reshape(dat_t,64,length(dat_t)/64);
data_offset=data+int16(repmat([(0:63)*200]',1,size(data,2)));
plot(data_offset')
