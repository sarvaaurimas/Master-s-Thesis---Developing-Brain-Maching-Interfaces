%misc variables
samp_rate = 30000; % smapling rate of open ephys system
pixels_per_m1 = 350;
pos_sample_rate1 = 50;
bin_size1=ceil(pixels_per_m1*0.025); % - it is around 2.5 cm x 2.5 cm per bin expressed


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
indeks = [1 4 6 35];
ratemaps = zeros(113,132,length(indeks)); % Now its 113 and 132 as that what size envir bin map is
posmaps = zeros(113,132,length(indeks)); % Now its 113 and 132 as that what size envir bin map is
ntimebins = round(length(xy1(:, 1))/50)-2; %+1 to make the bins fit the data and -2 as first 2s are discarded
max_time = 240500 %Use only first 20 mins of training data 
ntimebins = round(max_time/50)-2; %+1 to make the bins fit the data and -2 as first 2s are discarded
obsSpk = zeros(length(indeks), ntimebins);
% Some magic delay formula to account for diff in time started recording

% and 
delay = (syncStart-oeStart)/30 - 1000*postim(find(trig,1)); %in ms


%% run through all cells 
for i = 1:length(indeks)
    clf
    %get spike coordinates for a cell
    cell=indeks(i);
    disp([i length(indeks)])
    indz = sptime(klusters==cell); %30khz sampling
    %tc = double(indz)/30000;
    %changed cellNr_spike_coordinates1 to csc1.. csc1 sampled at 50Hz to
    %match with position sampling
    csc1 = round(double(indz)*50/30000 - delay*50/1000);% - 6*4096*50/1000); %  also some files were added where they shouldn't !
    csc1 = csc1(csc1>100); % 2 sec offset
%     csc1 = csc1(csc1<89800); %for30min trial
    csc1 = csc1(csc1<(length(x)-100)); %
%     csc1 = csc1(csc1<59250); %for20min trial


    %% start plotting
    if csc1(1)==0; csc1(1)=1; end
    csc1(1);
    edges = linspace(0, 240500, round(length(xy1(:, 1))/50)+1); % Bin boundaries, here multiples of 50 as position is sampled at 50hz so each bin corresponds to 1s
    [spike_bins, edges] = histcounts(csc1, edges); % Returns binned spike counts (-2 due to 2 sec offset applied earlier)
    spike_bins = spike_bins(3:end); % Discard first two seconds of time offset. 
    obsSpk(i, :) = spike_bins;
    disp('getting pos map');
    pos_map=Location_Map(xy1(1:max_time, :), bin_size1);
    disp('Getting spike map');
    sp_map=Spike_Map(xy1(1:max_time, :), csc1(find(csc1<=max_time)), bin_size1);
    %minx=(min(xy1(:, 1))+0.001);    maxx=max(xy1(:, 1));    miny=(min(xy1(:, 2))+0.001);    maxy=max(xy1(:, 2));
    %getting unsmoothened rate map for generating SAC and Fourier  2D spectrogram
    %sp_rmap=sp_map(ceil(miny/bin_size1):ceil(maxy/bin_size1), ceil((minx)/bin_size1):ceil(maxx/bin_size1))./pos_map(ceil(miny/bin_size1):ceil(maxy/bin_size1), ceil(minx/bin_size1):ceil(maxx/bin_size1));
    %Generating smoothened rate map using the adaptive smoothing:
    [smooth_pos_map, smooth_r_map]=apply_adaptive_smoothing(pos_map, sp_map, pos_sample_rate1);
    smooth_r_map = flipud(smooth_r_map); %Mirror flip it 180 degrees
    ratemaps(:, :, i) = smooth_r_map;
    posmaps(:, :, i) = smooth_pos_map;
    
end

nbins = round(max_time/50) - 2 % We disregarded first 2s
post = decode_calcBayesPost(obsSpk_2, ratemaps, posmaps(:, :, 1), 1); % Compute posteriors for each step 
max_post = reshape(max(max(post)), [nbins, 1]); % Calculate maximum posteriors in each step
row = zeros(length(max_post), 1);
col = zeros(length(max_post), 1);

for i = 1:length(max_post)
    curr_post = post(:, :, i);
    [row(i), col(i)] = find(post(:, :, i) == max_post(i)); % Get indexes of max posteriors in each step, sometimes throws errors if using too little cells
end
x = (row-0.5)*bin_size1;
y = (col-0.5)*bin_size1; % Predicted centers of bins
nbins = max_time/50 - 2 % We disregarded first 2s
x_real = mean(reshape(xy1(101:max_time, 1), [50, nbins]), 1);
y_real = mean(reshape(xy1(101:max_time, 2), [50, nbins]), 1);



