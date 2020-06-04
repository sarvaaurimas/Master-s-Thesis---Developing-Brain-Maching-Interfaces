
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

%% Compute x and y grids 
n_bins_x = ceil(max(xy1(:, 1)/bin_size1)); % number of bins in x axis
n_bins_y = ceil(max(xy1(:, 2)/bin_size1)); % number of bins in y axis
n_bins_z = 30;
bin_x_coord = ([1 : n_bins_x] -0.5)*bin_size1/(pixels_per_m1/100); % Coodinates of each bin
bin_y_coord = ([1 : n_bins_y] -0.5)*bin_size1/(pixels_per_m1/100);
[X_bin, Y_bin] = meshgrid(bin_x_coord, bin_y_coord); % A meshgrid for each bin center


disp('Encoding phase');
encoding_time = 7*60; % 5 minutes
n_tetrodes = 16;
n_electrodes = 1;
interspike_interval = 0.35; % 0.35 ms interspike interval
enc_raw_data = load_spikes(encoding_time); % The time in seconds is passed as and argument
[encoding_vectors, n_tetrode_spikes, sigmas] = compute_encoding_vectors_v2(enc_raw_data, n_tetrodes, n_electrodes, encoding_time, xy1, pixels_per_m1, interspike_interval);

training_pi_time = encoding_time;
B_x = [5, 5]; % the bandwidth for the X kde (6cm)
bw_volt = 5; % 24 mikroV for a start
B_dec = horzcat(ones(1, n_electrodes)*bw_volt, B_x);

pi_x = occupancy_x(X_bin, Y_bin, B_x, pixels_per_m1, training_pi_time, n_bins_y, n_bins_x, xy1);

tetrode_mu = n_tetrode_spikes/encoding_time; % Array of inidividual tetrode mu's
tetrode_lamb_x = tetrode_lambda_x(X_bin, Y_bin, B_x, n_bins_y, n_bins_x, encoding_vectors, n_tetrode_spikes, pi_x, n_tetrodes, tetrode_mu);
[tetrode_lamb_ax, X_bin_3d, Y_bin_3d, Z_bin_3d, tetrode_z_coord] = compute_tetrode_lambda_ax_v2(B_dec, encoding_vectors, tetrode_mu, pixels_per_m1, n_bins_y, n_bins_x, n_bins_z, pi_x, bin_size1, n_tetrodes, n_tetrode_spikes);


disp('Decoding phase');
decoding_time = 9*60; % 1 minute
dec_raw_data = enc_raw_data
%dec_filtered_data = filtered_fft(dec_raw_data, cutoff_freq); 

dt = 1; % 1s window

[log_post, decoded_position, n_dec_spikes, log_sum_lamb_ax_arr, dec_spikes] = compute_unsorted_posterior_position_precomp_v3(...
                                                           dec_raw_data, dt, n_bins_x, n_bins_y, n_tetrodes, ...
                                                       n_electrodes, pixels_per_m1, tetrode_lamb_x, bin_size1, ...
                                                       tetrode_z_coord, tetrode_lamb_ax, ...
                                                        interspike_interval, sigmas, X_bin, Y_bin, n_step_vel);
                                                    
% [log_post, decoded_position, n_dec_spikes, log_sum_lamb_ax_arr, dec_spikes] = compute_unsorted_posterior_position_precomp_v2(...
%                                                            dec_raw_data, dt, n_bins_x, n_bins_y, n_tetrodes, ...
%                                                        n_electrodes, pixels_per_m1, tetrode_lamb_x, bin_size1, ...
%                                                        tetrode_z_coord, tetrode_lamb_ax, ...
%                                                         interspike_interval, sigmas);
n_test_bins = size(decoded_position, 1);                                                               
x_real = (mean(reshape(xy1(1:floor(decoding_time/(dt))*(dt*50), 1), [dt*50, n_test_bins]), 1)/pixels_per_m1*100)';
y_real = (mean(reshape(xy1(1:floor(decoding_time/(dt))*(dt*50), 2), [dt*50, n_test_bins]), 1)/pixels_per_m1*100)';


[rmse, err] = xp_xr(decoded_position(:, 1), decoded_position(:, 2), x_real, y_real, size(decoded_position, 1));


