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


indeks = [1]; % Which cells to be plotted
ratemaps = zeros(113,132,size(indeks)); % Now its 113 and 132 as that what size envir bin map is

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
    tc = double(indz)/30000;
    %changed cellNr_spike_coordinates1 to csc1.. csc1 sampled at 50Hz to
    %match with position sampling
    csc1 = round(double(indz)*50/30000 - delay*50/1000);% - 6*4096*50/1000); %  also some files were added where they shouldn't !
    csc1 = csc1(csc1>100); % 2 sec offset
%     csc1 = csc1(csc1<89800); %for30min trial
    csc1 = csc1(csc1<(length(x)-100)); %
%     csc1 = csc1(csc1<59250); %for20min trial


    %%start plotting
    if csc1(1)==0; csc1(1)=1; end
    pos_map=Location_Map(xy1, bin_size1);
    sp_map=Spike_Map(xy1, csc1, bin_size1);
   
    minx=(min(xy1(:, 1))+0.001);    maxx=max(xy1(:, 1));    miny=(min(xy1(:, 2))+0.001);    maxy=max(xy1(:, 2));
    %getting unsmoothened rate map for generating SAC and Fourier  2D spectrogram
    sp_rmap=sp_map(ceil(miny/bin_size1):ceil(maxy/bin_size1), ceil((minx)/bin_size1):ceil(maxx/bin_size1))./pos_map(ceil(miny/bin_size1):ceil(maxy/bin_size1), ceil(minx/bin_size1):ceil(maxx/bin_size1));
    %Generating smoothened rate map using the adaptive smoothing:
    %[smooth_pos_map, smooth_r_map]=apply_adaptive_smoothing(pos_map, sp_map, pos_sample_rate1);
    %smooth_r_map = flipud(smooth_r_map); %Mirror flip it 180 degrees
    Draw_RMap(sp_rmap);
end
    
