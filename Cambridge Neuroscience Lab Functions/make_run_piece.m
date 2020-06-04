clf
disp(['cell number ' num2str(cell)])
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