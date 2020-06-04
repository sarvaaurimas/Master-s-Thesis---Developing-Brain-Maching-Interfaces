
%misc variables
samp_rate = 30000; % smapling rate of open ephys system
pixels_per_m1 = 350;
pos_sample_rate1 = 50;
bin_size1=ceil(pixels_per_m1*0.025); % - it is around 2.5 cm x 2.5 cm per bin expressed
tBin = 0.25; %0.25s time bin for spikes
cutoff_freq = 300; % To filter all LFPs



%Load all the vars
trial_folder = 'C:\Users\Aurimas\Desktop\Uniko stuff\4th year\4th year project\files';
posfile = fullfile(trial_folder,'b_20180823_110359.bin');
[xys,postim,trig] = clean_positions(posfile); %loads xy coordinates and timestamp for respective coordinates'
sptime = readNPY('spike_times.npy');
%klusters = readNPY('spike_clusters.npy');
%pc_features = readNPY('pc_features.npy');
%pc_features_ind = readNPY('pc_feature_ind.npy');

sync_mess = fullfile(trial_folder,'sync_messages.txt');
ch_tmpstmp = fullfile(trial_folder,'timestamps_ttl.npy');
datatims = fullfile(trial_folder,'timestamps_cont.npy');



%syncTimer has timestamps where trig pulse was read
syncTimer=double(readNPY(ch_tmpstmp)); syncStart=syncTimer(1); syncTimer=(syncTimer-syncTimer(1))/30000;
%oeTimer contains timestamps for each ephys sample
oeTimer=double(readNPY(datatims)); oeTimer=(oeTimer-oeTimer(1))/30000;
oeStart = clean_startime(sync_mess); %start of sample clock of oe system


% Calculate coordinates
a1=xys(:,1);
a2=xys(:,2);
b1=xys(:,4);
b2=xys(:,5);
x=0.5*(a1+b1); y=0.5*(a2+b2);
angles=atan2(a1-b1,a2-b2)*180/pi;
xy1 = [x,y];


%indeks = [1 4 6 35 69 91 94 107 117 142 163 210 214 224 232 244 245]; % Which cells to be plotted
%upd_indeks = [2, 3, 5, 7:34, 36:68, 70:90, 92, 93, 95:106, 108:116,
%118:141, 143:162, 164:209, 211, 212, 213, 215:223, 225:231, 233:243,
%246:256]
%indeks = [1:256]
%indeks = [1]
%indeks = [1 4 6 35 69 91 94 107 117 ];
%indeks = [1]
%max_time_test = 240500/4 %Use only first 20 mins of training data 
%max_time = 240500/4 %Use all 80 mins, this is in 50 Hz sampled
max_time = 240500;
%ntimebins = round(max_time/50)-3; % - 3 for div/4, -2 for full
%obsSpk = zeros(length(indeks), ntimebins);
% Some magic delay formula to account for diff in time started recording
delay = (syncStart-oeStart)/30 - 1000*postim(find(trig,1)); %in ms

%changed cellNr_spike_coordinates1 to csc1.. csc1 sampled at 50Hz to
%match with position sampling
% Line below for just cell 1
csc1 = round(double(sptime)*50/30000 - delay*50/1000);% - 6*4096*50/1000); %  also some files were added where they shouldn't !
%csc1 = round(double(sptime(klusters==1))*50/30000 - delay*50/1000);% - 6*4096*50/1000); %  also some files were added where they shouldn't !
%csc1 = round(double(sptime)*50/30000 - delay*50/1000);% - 6*4096*50/1000); %  also some files were added where they shouldn't !
csc1 = csc1(csc1>100); % 2 sec offset
%csc1 = csc1(csc1<89800); %for30min trial
csc1 = csc1(csc1<(length(x)-100)); %
%csc1 = csc1(csc1<59250); %for20min trial

if size(csc1, 1)== 0
% Do nothing
else
    %% start plotting
    if csc1(1)==0; csc1(1)=1; end
end

n_bins_x = ceil(max(xy1(:, 1)/bin_size1)); % number of bins in x axis
n_bins_y = ceil(max(xy1(:, 2)/bin_size1)); % number of bins in y axis
bin_x_coord = ([1 : n_bins_x] -0.5)*bin_size1/(pixels_per_m1/100); % Coodinates of each bin
bin_y_coord = ([1 : n_bins_y] -0.5)*bin_size1/(pixels_per_m1/100);
[X_bin, Y_bin] = meshgrid(bin_x_coord, bin_y_coord); % A meshgrid for each bin center

disp('Encoding phase');
encoding_time = 3*60; % 5 minutes
n_tetrodes = 8;
n_electrodes = 1;
spike_distance = 7;
enc_raw_data = load_spikes(encoding_time); % The time in seconds is passed as and argument
[encoding_vectors, n_tetrode_spikes] = compute_encoding_vectors(enc_raw_data, n_tetrodes, n_electrodes, encoding_time, xy1, spike_distance, pixels_per_m1);

disp('Decoding phase');
decoding_time = 0.3*60; % 1 minute
dec_raw_data = load_spikes(decoding_time);
%dec_filtered_data = filtered_fft(dec_raw_data, cutoff_freq); 


dt = 1; % 1s window
B_x = [6, 6]; % the bandwidth for the X kde (6cm)
spike_coords = xy1(csc1(csc1 < max_time), :)/pixels_per_m1*100; % all true spike coords 
N = size(csc1(csc1 < max_time), 1); % Number of spikes fired
T = max_time/50; % Total training time
mu = N / T; % Average firing rate per second
[lamb_x, pi_x] = lambda_x(X_bin, Y_bin, B_x, spike_coords, mu, pixels_per_m1, max_time, n_bins_y, n_bins_x, xy1);
bw_volt = 24; % 24 mikroV for a start
B_dec = horzcat(ones(1, n_electrodes)*bw_volt, B_x);

[log_post, decoded_position, n_enc_spikes, log_sum_lamb_ax_arr] = compute_unsorted_posterior_position(dec_raw_data, dt, encoding_vectors, pi_x,...
                                                                   n_bins_x, n_bins_y, n_tetrodes, spike_distance, ...
                                                                   n_electrodes, X_bin, Y_bin, B_dec, mu, xy1, ...
                                                                   pixels_per_m1, max_time, n_tetrode_spikes, lamb_x, bin_size1);

n_test_bins = size(decoded_position, 1);                                                               
x_real = (mean(reshape(xy1(1:floor(decoding_time/(dt))*(dt*50), 1), [dt*50, n_test_bins]), 1)/pixels_per_m1*100)';
y_real = (mean(reshape(xy1(1:floor(decoding_time/(dt))*(dt*50), 2), [dt*50, n_test_bins]), 1)/pixels_per_m1*100)';


[rmse, err] = xp_xr(decoded_position(:, 1), decoded_position(:, 2), x_real, y_real, size(decoded_position, 1));

n_slice = 3;
%surf(reshape(X_bin(n_slice, :, :), size(X_bin, 2),size(X_bin, 3)) , reshape(Z_bin(3, :, :), size(X_bin, 2),size(X_bin, 3)), reshape(lamb_ax(3, :, :), size(X_bin, 2),size(X_bin, 3)));


