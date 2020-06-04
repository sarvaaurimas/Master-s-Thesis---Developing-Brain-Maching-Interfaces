%misc variables
samp_rate = 30000; % smapling rate of open ephys system
pixels_per_m1 = 350;
pos_sample_rate1 = 50;
bin_size1=ceil(pixels_per_m1*0.025); % - it is around 2.5 cm x 2.5 cm per bin expressed
tBin = 1; %1s time bin for spikes

%Load all the vars
trial_folder = 'C:\Users\Aurimas\Desktop\Uniko stuff\4th year\4th year project\files';
posfile = fullfile(trial_folder,'b_20180823_110359.bin');
[xys,postim,trig] = clean_positions(posfile); %loads xy coordinates and timestamp for respective coordinates'
sptime = readNPY('spike_times.npy');
klusters = readNPY('spike_clusters.npy');

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
indeks = [1:256]
%indeks = [1, 4, 6]
%indeks = [12];
%indeks = [1]
ratemaps = zeros(113,132,0); % Now its 113 and 132 as that what size envir bin map is, 113x130 if /4
max_time_test = 240500/4 %Use only first 20 mins of training data 
max_time = 240500 %Use all 80 mins, this is in 50 Hz sampled
ntimebins = round(max_time/50)-2; % - 3 for div/4, -2 for full
obsSpk = zeros(length(indeks), ntimebins);
% Some magic delay formula to account for diff in time started recording
delay = (syncStart-oeStart)/30 - 1000*postim(find(trig,1)); %in ms


%% run through all cells 
%Calculate Position map once
disp('getting pos map');
pos_map=Location_Map(xy1(1:max_time, :), bin_size1);
for i = 1:length(indeks)
    clf
    %get spike coordinates for a cell
    cell=indeks(i);
    disp([i length(indeks)]);
    indz = sptime(klusters==cell); %30khz sampling, times when that cell spiked
    %tc = double(indz)/30000;
    %changed cellNr_spike_coordinates1 to csc1.. csc1 sampled at 50Hz to
    %match with position sampling
    csc_unsorted = round(double(indz)*50/30000 - delay*50/1000);% - 6*4096*50/1000); %  also some files were added where they shouldn't !
    csc_unsorted = csc_unsorted(csc_unsorted>100); % 2 sec offset
%     csc1 = csc1(csc1<89800); %for30min trial
    csc_unsorted = csc_unsorted(csc_unsorted<(length(x)-100)); %
%     csc1 = csc1(csc1<59250); %for20min trial

    if size(csc_unsorted, 1)== 0
    % Do nothing
    else 
    %% start plotting
    if csc_unsorted(1)==0; csc_unsorted(1)=1; end
    csc_unsorted(1);

    
    [spkCnt, winLim] = get_SpkCnt(csc_unsorted, [1, max_time], 50); %Try to get spike counts with new function
    obsSpk(i, :) = spkCnt(3:end);
    disp('Getting spike map');
    sp_map=Spike_Map(xy1(1:max_time, :), csc_unsorted(find(csc_unsorted<=max_time)), bin_size1); % in each bin is the count of spikes in that bin
   
    %Generating smoothened rate map using the adaptive smoothing:
    disp('smoothing');
    [smooth_pos_map, smooth_r_map]=apply_adaptive_smoothing(pos_map, sp_map, pos_sample_rate1);
    [norm_pos_map, norm_r_map] = normalize_maps(pos_map, sp_map, pos_sample_rate1);
    %smooth_r_map = flipud(smooth_r_map); %Mirror flip it 180 degrees
    ratemaps(:, :, i) = smooth_r_map;
    ratemaps = cat(3, ratemaps, smooth_r_map);
    end 
end
posmap= smooth_pos_map;
n_cells = size(ratemaps, 3);
n_test_bins = round(max_time_test/50) - 3;

[x, y] = decode_posPred(n_test_bins, obsSpk, ratemaps, posmap, tBin, bin_size1, pixels_per_m1, n_cells);
x_real = (mean(reshape(xy1(101:floor(max_time_test/50)*50, 1), [50, n_test_bins]), 1)/pixels_per_m1*100)';
y_real = (mean(reshape(xy1(101:floor(max_time_test/50)*50, 2), [50, n_test_bins]), 1)/pixels_per_m1*100)';

figure(1)
err = xp_xr(x, y, x_real, y_real, 1000);
%figure(2)
%plot_rmse(obsSpk_2(:, 1:n_test_bins), ratemaps, posmaps, x_real, y_real, 1200, bin_size1, pixels_per_m1); 

%% FUNCTION GUIDE
% decode_posPred - decodes position from data
% 
% PLOTTING FUNCTIONS
% xp_xr - plots real vs predicted and rms
% corellogram - plots a correlogram between xp and xr
% plot_vConstVar - plots how rms varies with different K constant
% plot_cong_err - plots scatter plot and best fit line of confidence vs rms
% err_big_data - computes rms for data too big to be analysed in one go
% plot_posterior - plots the posterior function distribution @ the actual
% current position


% scatter(max_post(1:n_test_bins/5),err(1:n_test_bins/5))
% ylabel('RMSE, cm')
% xlabel('Confidence, max P(x|n)')

%rmse = plot_vConstVar(n_test_bins, obsSpk_2, ratemaps, posmaps, x_real, y_real, bin_size1, pixels_per_m1, tBin);

%figure(1)
%corellogram(x, x_real);
% nbins = round(max_time/50) - 2 % We disregarded first 2s
% post = decode_calcBayesPost(obsSpk_2, ratemaps, posmaps(:, :, 1), 1); % Compute posteriors for each step 
% max_post = reshape(max(max(post)), [nbins, 1]); % Calculate maximum posteriors in each step
% row = zeros(length(max_post), 1);
% col = zeros(length(max_post), 1);
% 
% for i = 1:length(max_post)
%     curr_post = post(:, :, i);
%     [row(i), col(i)] = find(post(:, :, i) == max_post(i)); % Get indexes of max posteriors in each step, sometimes throws errors if using too little cells
% end
% y = (row-0.5)*bin_size1/pixels_per_m1*100;
% x = (col-0.5)*bin_size1/pixels_per_m1*100; % Predicted centers of bins in pixels
% x_real = (mean(reshape(xy1(101:floor(max_time/50)*50, 1), [50, nbins]), 1)/pixels_per_m1*100)';
% y_real = (mean(reshape(xy1(101:floor(max_time/50)*50, 2), [50, nbins]), 1)/pixels_per_m1*100)';
% figure(1)
% corellogram(x, x_real);
% figure(2)
% err = xp_xr(x, y, x_real, y_real, 1000);
% 
% plot_rmse(obsSpk_2, ratemaps, posmaps(:, :, 1), x_real, y_real, 1200, bin_size1, pixels_per_m1); 



