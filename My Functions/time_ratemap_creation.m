function [timestamps, elapsed_time] = time_ratemap_creation(max_time)
n = 10;
cells =5;
timestamps = linspace(100, max_time, n);
elapsed_time = zeros(1,n);

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


indeks = [1:cells]

ratemaps = zeros(113,132,0); % Now its 113 and 132 as that what size envir bin map is, 113x130 if /4
max_time_test = 240500/4 %Use only first 20 mins of training data 
max_time = 240500 %Use all 80 mins, this is in 50 Hz sampled
ntimebins = round(max_time/50)-2; % - 3 for div/4, -2 for full
obsSpk = zeros(length(indeks), ntimebins);
% Some magic delay formula to account for diff in time started recording
delay = (syncStart-oeStart)/30 - 1000*postim(find(trig,1)); %in ms


%% run through all cells 
%Calculate Position map once




for h = 1:n
    disp('Number of main loop iterations');
    disp(h);
    max_time = timestamps(h);
    disp('getting pos map');
    pos_map=Location_Map(xy1(1:max_time, :), bin_size1);
    tic
    
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

    
    disp('Getting spike map');
    sp_map=Spike_Map(xy1(1:max_time, :), csc_unsorted(find(csc_unsorted<=max_time)), bin_size1); % in each bin is the count of spikes in that bin
   
    %Generating smoothened rate map using the adaptive smoothing:
    disp('smoothing')
    [smooth_pos_map, smooth_r_map]=apply_adaptive_smoothing(pos_map, sp_map, pos_sample_rate1);
    end
    
    end
    
    elapsed_time(h) = toc;
    disp('Elapsed time')
    disp(elapsed_time(h));
end
elapsed_time = elapsed_time/cells;
end
