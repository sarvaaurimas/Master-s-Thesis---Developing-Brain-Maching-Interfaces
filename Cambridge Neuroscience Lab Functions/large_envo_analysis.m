%% this file contains workflow required to process open ephys data in large
%% (and small) environmnets. General flow is like this:
% Ephys recordings starts
% Tracking starts
% Tracking finishes
% Ephys recordings finishes
% Transfer data to storage folder
% Filter binary data
% Sort binary data
% Curate binary data
% Extract position data
% Sync both of them
% Create ratemaps
% Output all of it to matlab (python) format and prepare for further analysis

%misc variables
samp_rate = 30000; % smapling rate of open ephys system

%% Filtering binary data
%move position data for a trial into main recording folder for the trial
%open the main recording folder for this particular trial
trial_folder = 'D:\data\R3\a_2018-05-28_13-27-47';
posfile = fullfile(trial_folder,'a_20180528_132704.bin');

sync_mess = fullfile(trial_folder,'experiment1\recording1\sync_messages.txt');
ch_states = fullfile(trial_folder,'experiment1\recording1\events\Rhythm_FPGA-100.0\TTL_1\channel_states.npy');
ch_tmpstmp = fullfile(trial_folder,'experiment1\recording1\events\Rhythm_FPGA-100.0\TTL_1\timestamps.npy');
databin = fullfile(trial_folder,'experiment1\recording1\continuous\Rhythm_FPGA-100.0\continuous.dat');
datatims = fullfile(trial_folder,'experiment1\recording1\continuous\Rhythm_FPGA-100.0\timestamps.npy');
datafile = fullfile(trial_folder,'spikes.bin');



[xys,postim,trig] = clean_positions(posfile); %loads xy coordinates and timestamp for respective coordinates
% % % plot(xys(:,1),xys(:,2),'.')
% % % [px, py] = ginput(4);
% % % xxXx=inpolygon(xys(:,1),xys(:,2),px,py);
% % % xys(xxXx==0,:)=NaN;
% need to interpolate to 50hz here





cd(trial_folder);

Hd = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',500,'HalfPowerFrequency2',6000,'SampleRate',30000);

fib = fopen(databin); fib2 = fopen(datafile,'w');% fib3 = fopen(fullfile(trial_folder,'lfp.bin'),'w'); %open files
fll = dir(databin);
numbz = floor(fll.bytes/2/70/1e7);
for ji = 1:numbz
datt = fread(fib,[70 1e7],'*int16','ieee-le'); %64 neural signals + 6 aux inputs
datt=0.195*double(datt(1:64,:));
datt = bsxfun(@minus,datt,mean((datt),2));
datt(1:32,:) = bsxfun(@minus,datt(1:32,:),median(datt(1:32)));
datt(33:64,:) = bsxfun(@minus,datt(33:64,:),median(datt(33:65)));
for i = 1:64
    datt(i,:) = int16(filtfilt(Hd,datt(i,:)));    
end
fwrite(fib2,datt,'int16');
end
datt = fread(fib,[70 Inf],'*int16','ieee-le'); %64 neural signals + 6 aux inputs
datt=0.195*double(datt(1:64,:));
datt = bsxfun(@minus,datt,mean(datt,2));
datt(1:32,:) = bsxfun(@minus,datt(1:32,:),median(datt(1:32)));
datt(33:64,:) = bsxfun(@minus,datt(33:64,:),median(datt(33:65)));
for i = 1:64
    datt(i,:) = filtfilt(Hd,datt(i,:));    
end
fwrite(fib2,datt,'int16');


% [spikes lfp] = clean_data(datt); %spikes - 30kHz. LFP - 1kHz
 fclose(fib2); %fwrite(fib3,lfp); fclose(fib3); %close files

% xxyy = xys(:,1:2); save('pozicija.mat','xxyy');
syncTimer=double(readNPY(ch_tmpstmp)); syncStart=syncTimer(1); syncTimer=(syncTimer-syncTimer(1))/30000;
oeTimer=double(readNPY(datatims)); oeTimer=(oeTimer-oeTimer(1))/30000;

oeStart = clean_startime(sync_mess); %sample clock of oe system

%% KILOSORT
masterstva(trial_folder,datafile);

% REMOVE OUTSIDE STUFF%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plotting
sptime = readNPY('spike_times.npy');
klusters = readNPY('spike_clusters.npy');
figFol = fullfile(trial_folder,'figures');
if ~exist(figFol)
    mkdir(figFol)
end
cd(figFol);
a1=xys(:,1);
a2=xys(:,2);
b1=xys(:,4);
b2=xys(:,5);
x=a1; y=a2;
angles=atan2(a1-b1,a2-b2)*180/pi;
cut1 = klusters;
output1 = sptime;


klus = unique(klusters);
indeks=[];
for i = 1:length(klus)
    ind = sptime(klusters==klus(i));
    if length(ind)>100 %remove separate spikes
        indeks = [indeks klus(i)];
    end
end

%% preps
pixels_per_m1 = 350;
pos_sample_rate1 = 50;
bin_size1=ceil(pixels_per_m1*0.025); % - it is around 2.5 cm x 2.5 cm per bin expressed
bin_size1=ceil(pixels_per_m1*0.03); % - 5x5 coz bad sampling

xy1 = [x,y];
spbns=0:45; %range where to check speed (in cm/s)
bb = (1/5)*ones(1,5);
aa=1;
xxx=filter(bb,aa,x);
yyy=filter(bb,aa,y);
speed = [0; sqrt( diff(xxx).^2 + diff(yyy).^2 )];
speed(speed>6.5)=0;
spy=speed*50/2.3;%50hz sampling, 2.3px per cm
gap=50;apo_cr = repmat((0:gap:(64-1)*gap)',1,151);
% gap=50;apo_cr = repmat((0:gap:(192-1)*gap)',1,151);
[count, coo] = hist(angles,72);
angles(angles==0)=NaN;
[an_sp, Centers] = hist(spy,spbns);


delay=0;
% first - when recording started
%second - when ttracking started
delay = 11530; %in miliseconds
delay = (syncStart-oeStart)/30 - 1000*postim(find(trig,1)); %in ms
% delay = delay*32; %32kHz sampling rate 
% indeks = [58,64,191,270,33,116,16,30,118,42,79,136,159,229,106,140,62,141,105,18,240,196,88,18,9,15,175,91,24,180,83,89,193,106,62,207,138,58,86,51,229,44,276,68,211,274,63,219,224,172,47,64];
% indeks = [208 72 85 81 93 103 88 29 40 83 210 255 25 57 113 190 6 45];
% indeks = [23 100 131 235 123 79 108 1 140 219 0 18 89 107];
%% run through all cells
for i = 1:length(indeks)
    clf
    %get spike coordinates for a cell
    cell=indeks(i);
    disp([i length(indeks)])
    indz = sptime(klusters==cell); %32khz sampling
    tc = double(indz)/30000;
    %changed cellNr_spike_coordinates1 to csc1.. csc1 sampled at 50Hz to
    %match with position sampling
    csc1 = round(double(indz)*50/30000 - delay*50/1000);% - 6*4096*50/1000); %  also some files were added where they shouldn't !
    csc1 = csc1(csc1>100); % 2 sec offset
%     csc1 = csc1(csc1<89800); %for30min trial
    csc1 = csc1(csc1<(length(x)-100)); %
%     csc1 = csc1(csc1<59250); %for20min trial
    %% start plotting
    if csc1(1)==0; csc1(1)=1; end
    
    subplot(5,5,[1 2 6 7]);
    %plot spatial fields here

    
    pos_map=Location_Map(xy1, bin_size1);
    sp_map=Spike_Map(xy1, csc1, bin_size1);
    
    minx=(min(xy1(:, 1))+0.001);    maxx=max(xy1(:, 1));    miny=(min(xy1(:, 2))+0.001);    maxy=max(xy1(:, 2));
    %getting unsmoothened rate map for generating SAC and Fourier  2D spectrogram
    sp_rmap=sp_map(ceil(miny/bin_size1):ceil(maxy/bin_size1), ceil((minx)/bin_size1):ceil(maxx/bin_size1))./pos_map(ceil(miny/bin_size1):ceil(maxy/bin_size1), ceil(minx/bin_size1):ceil(maxx/bin_size1));
    %Generating smoothened rate map using the adaptive smoothing:
    [smooth_pos_map, smooth_r_map]=apply_adaptive_smoothing(pos_map, sp_map, pos_sample_rate1);
    Draw_RMap(smooth_r_map);
    axis equal;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    %text(2, 60, [ num2str(roundn(max(max(SpRmap(isnan(SpRmap)~=1))), -1)) ' Hz']);
    set(gca, 'visible', 'off');
    set(gca,'YDir','normal')

   
    
    
    subplot(5,5,[11 12 16 17]);
    %plot spike map here
    plot(x,y,'.k')
    hold on
    plot(x(ceil(csc1)),y(ceil(csc1)),'r.')
    axis equal;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    %text(2, 60, [ num2str(roundn(max(max(SpRmap(isnan(SpRmap)~=1))), -1)) ' Hz']);
    set(gca, 'visible', 'off');
    
    subplot(5,5,4);
    %plot SAC here
    sacSmooth1=get_smooth_SAC(sp_rmap);
    Draw_RMap(sacSmooth1);
    axis equal;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    %text(2, 60, [ num2str(roundn(max(max(SpRmap(isnan(SpRmap)~=1))), -1)) ' Hz']);
    set(gca, 'visible', 'off');
    title('Spatial AC')
    
    subplot(5,5,[3 8 13 18])
    %plot waveforms here
%     sumarum = -squeeze(oot.waveFormsMean(i,:,:));
%     plot((-sumarum(:,25:65)+apo_cr(:,40:80))','k')
% % % %     ylim([-300 19000])
%      xlim([0 40])
%      delete(findall(ancestor(p,'figure'),'HandleVisibility','off','type','line','-or','type','text'));
sumarum = apo_cr(:,40:80);     

     %get WF parameters
     try
     chn = 1:64;
     amps = sumarum(:,41);
     f = fit(chn',amps,'gauss1');
     rr=coeffvalues(f);
     
     amplitude = max(max(sumarum));
     center = rr(2);
     tetrode_number = ceil(center/4);
     depth = (center * 10)-10;     
     %peak-troph
     cno = round(center);
     if cno<1; cno=1; end
     if cno>64; cno=64; end
     chaa = sumarum(cno,:);     
     fwhm = sum(sumarum(:,41)>0.3*amplitude);     
     peatrop = (find(chaa==min(chaa))-find(chaa==max(chaa)))*32.25;
     catch
        rr = [0 0 0]; 
     amplitude = max(max(sumarum));
     center = find(sumarum(:,41) == max(sumarum(:,41)));
     depth = (center * 10)-10;
     fwhm = sum(sumarum(:,41)>0.5*amplitude);
     %peak-troph
     cno = round(center);
     if cno<1; cno=1; end
     if cno>64; cno=64; end
     chaa = sumarum(cno,:);
     peatrop = (find(chaa==min(chaa))-find(chaa==max(chaa)))*30;
     end
     
     amp_val = [];
     centras = find(sumarum(:,41) == max(sumarum(:,41)));
     
     
    subplot(5,5,[21:25])
    %plot firing dependence on time here
% 	delete(findall(ancestor(p,'figure'),'HandleVisibility','off','type','text'));
%     plot(inds,amp_val,'k.')
plot(csc1,rand(1,length(csc1)),'k.')
    xlim([0 length(x)])
    subplot(5,5,9)
    %plot polar plot here
    [count1, coo1] = hist(angles(csc1),72);
    p=polar((coo)/57.3,smooth(count1./count)','r');
    title('Polar plot')
    delete(findall(ancestor(p,'figure'),'HandleVisibility','off','type','text'));
     axis equal;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    %text(2, 60, [ num2str(roundn(max(max(SpRmap(isnan(SpRmap)~=1))), -1)) ' Hz']);
    set(gca, 'visible', 'off');
%     
    
    if length(tc) > 10000
        tc=tc(1:10000);
    end
        subplot(5,5,5)
        %plot autocorr here with windows of 10
        [bins cou] = interspike_histogram(tc, tc, 10 );
        bar(bins,cou)
        title([num2str(sum(cou(52:60))) '/' num2str(length(tc))])
        xlim([-11 11])
        subplot(5,5,10)
        %plot autocorr here with windows of 500
        [bins cou] = interspike_histogram(tc, tc, 500 );
        bar(bins,cou)
        xlim([-550 550])
   
    subplot(5,5,[14 15])
    %plot speed profile here
    
    cell_speed = spy(csc1);
        [cell_sp, Centers] = hist(cell_speed,spbns);
        norm_sp = cell_sp./an_sp;
        bar(Centers,norm_sp);
        xlim([0 40])
        ylim([0 max(norm_sp(1:30))*1.2])
        title('Rate vs Speed')
        
        subplot(5,5,19)
        set(gca, 'visible', 'off');
        text(0,-0.2,'Number of spikes')
        text(0,0,'Mean firing rate')
        text(0,0.2,'Peak-troph (us)')
        text(0,0.4,'Amplitude')
        text(0,0.6,'tetrode number')
        text(0,0.8,'Cell number')
        text(0,1,'Trial')
        
        subplot(5,5,20)
        set(gca, 'visible', 'off');
        text(0,-0.2,num2str(length(tc)))        
        text(0,0,num2str(length(csc1)/1800))
        text(0,0.2,num2str(peatrop))
        text(0,0.4,num2str(amplitude))
        text(0,0.6,num2str(tetrode_number))
        text(0,0.8,num2str(cell))
        text(0,1,posfile)
    
   
    %plot textual info here: firing rate, max and mean firing rates, amplitude,
    %peak-troph timing, position of the peak and some other parameters, cell number
    %also need plot for directional info, 
    %do i need trajectory+dots as well?
    
    
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 0.9 0.9]);
%     guvna(indeks(i),klusters,sptime,xx,yy);
%   if depth <1800 && length(tc) >150
    print(['cell' num2str(cell)],'-dpdf')
%   end
    clear rr sumarum cell_speed norm_sp bins cou coo1 count1
    
end














function [data lfp] = clean_data(datt) %filter data
Hd = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',500,'HalfPowerFrequency2',6000,'SampleRate',30000);
hpFilt = designfilt('highpassiir','FilterOrder',8,'PassbandFrequency',1,'PassbandRipple',0.2,'SampleRate',1000);
lpFilt = designfilt('lowpassiir','FilterOrder',8,'PassbandFrequency',150,'PassbandRipple',0.2,'SampleRate',1000);
data = 0.195*double(datt(1:64,:)); %not using accelerometer
lfp = downsample(data',30)';
for i = 1:64
    lfp(i,:) = filtfilt(lpFilt,double(lfp(i,:)));
    lfp(i,:) = filtfilt(hpFilt,double(lfp(i,:)));
    data(i,:) = filtfilt(Hd,data(i,:));    
end
data(1:32,:) = bsxfun(@minus,data(1:32,:),mean(data(1:32,:),2)); %remove offsets HS1
data(33:64,:) = bsxfun(@minus,data(33:64,:),mean(data(33:64,:),2)); %remove offsets HS2
data(1:32,:) = bsxfun(@minus,data(1:32,:),median(data(1:32,:))); %remove noise HS1
data(33:64,:) = bsxfun(@minus,data(33:64,:),median(data(33:64,:))); %remove noise HS2
data = int16(data);
end







