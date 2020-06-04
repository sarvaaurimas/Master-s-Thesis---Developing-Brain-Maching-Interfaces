function [upd_ratemaps, upd_obsSpk]= updateRatemaps(ratemaps, obs_spk, newIndexes, max_time)
% Function that takes in already existing ratemaps, new indexes of
% cells and returns updated ratemaps

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


%upd_ratemaps = zeros(113,130,length(newIndexes) + size(ratemaps, 3)); % Now its 113 and 132 as that what size envir bin map is, 113x130 if /4
ntimebins = round(max_time/50)-3; % - 3 for div/4, -2 for full
% Some magic delay formula to account for diff in time started recording
delay = (syncStart-oeStart)/30 - 1000*postim(find(trig,1)); %in ms
%upd_ratemaps(:, :, 1:size(ratemaps, 3)) = ratemaps;
%upd_obsSpk = zeros(length(newIndexes)+size(ratemaps, 3), ntimebins);
%upd_obsSpk(1:size(ratemaps, 3), :) = obs_spk;
upd_ratemaps = ratemaps;
upd_obsSpk = obs_spk;
%% run through all cells 
disp('getting pos map');
pos_map=Location_Map(xy1(1:max_time, :), bin_size1);
for i = 1:length(newIndexes)
    clf
    %get spike coordinates for a cell
    cell=newIndexes(i);
    disp([i length(newIndexes)])
    indz = sptime(klusters==cell); %30khz sampling
    %tc = double(indz)/30000;
    %changed cellNr_spike_coordinates1 to csc1.. csc1 sampled at 50Hz to
    %match with position sampling
    csc1 = round(double(indz)*50/30000 - delay*50/1000);% - 6*4096*50/1000); %  also some files were added where they shouldn't !
    csc1 = csc1(csc1>100); % 2 sec offset
%     csc1 = csc1(csc1<89800); %for30min trial
    csc1 = csc1(csc1<(length(x)-100)); %
%     csc1 = csc1(csc1<59250); %for20min trial

    if size(csc1, 1)== 0
    % Do nothing
    else 
    %% start plotting
    
    if csc1(1)==0; csc1(1)=1; end
    csc1(1);
    [spkCnt, winLim] = get_SpkCnt(csc1, [1, max_time], 50); %Try to get spike counts with new function
    upd_obsSpk = cat(1, upd_obsSpk, spkCnt(3:end)'); %Transpose to make it 
    %upd_obsSpk(i+size(ratemaps, 3), :) = spkCnt(3:end);
    
    disp('Getting spike map');
    sp_map=Spike_Map(xy1(1:max_time, :), csc1(find(csc1<=max_time)), bin_size1);
   
    %Generating smoothened rate map using the adaptive smoothing:
    disp('smoothing');
    [smooth_pos_map, smooth_r_map]=apply_adaptive_smoothing(pos_map, sp_map, pos_sample_rate1);
    %smooth_r_map = flipud(smooth_r_map); %Mirror flip it 180 degrees
    %upd_ratemaps(:, :, i+size(ratemaps, 3)) = smooth_r_map;
    upd_ratemaps = cat(3, upd_ratemaps, smooth_r_map);
    Draw_RMap(smooth_r_map);
    end
end



