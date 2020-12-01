function plot_HP_ArmReplay_triggered_spatialprecision(thisdir,label)

load(thisdir,[label '_replay_singlebothjoint'],[label '_replay_replayarm'], ...
   'hp_cells','hpinterneurons',[label '_replay_seqmean'],[label '_replay_seqtimes'],[label '_replay_stnd'],...
    [label '_Cand_sig_modu_include'],[label '_pSSDarm'],[label '_replay_postseqindex'],'spikedata','armposindex')


eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['seqmean = ' label '_replay_seqmean;'])
eval(['seqtimes = ' label '_replay_seqtimes;'])
eval(['stnd = ' label '_replay_stnd;'])
eval(['postseqindex = ' label '_replay_postseqindex;'])
eval(['p_SSD = ' label '_pSSDarm;'])

clear([label '_replay_singlebothjoint'],[label '_replay_stnd'],[label '_replay_replayarm'], ...
   [label '_replay_seqmean'],[label '_replay_seqtimes'],...
    [label '_Cand_sig_modu_include'],[label '_pSSDarm'],[label '_replay_postseqindex'])

hp = hp_cells(~ismember(hp_cells,hpinterneurons)); 
replayarm(singlebothjoint==3) = NaN;
clear singlebothjoint
newcol = [75 0 130;255 130 0;34 139 34]/255;
[~,seq] = histc(seqmean,1:size(armposindex,1));

spks = spikedata(ismember(spikedata(:,2),hp),:);
postseqindex = [postseqindex;length(seq)+1];

sp = zeros(size(armposindex,1),length(hp));
ra = zeros(size(armposindex,1),1);
toadd = [0 find(armposindex(:,1),1,'last') find(armposindex(:,2),1,'last')];
for ii = 1:length(replayarm)
    if ~isnan(replayarm(ii))
        repind = postseqindex(ii):(postseqindex(ii+1)-1);
        reppos = seq(repind)+toadd(replayarm(ii));
        for icell = 1:length(hp)
           spks2 = spks(ismember(spks(:,2),hp(icell)),:); 
           dat = spks2(spks2(:,1)>=stnd(ii,1) & spks2(:,1)<stnd(ii,2),:);
           if ~isempty(dat)       
                [~,c] = histc(dat(:,1),seqtimes(repind));
                h = histc(c,unique(c));
                sp(reppos(unique(c)),icell) = sp(reppos(unique(c)),icell)+h;                
           end
        end
        h = histc(reppos,unique(reppos));
        ra(unique(reppos),1) = ra(unique(reppos),1)+h;
    end
end

% Filter=fspecial('gaussian',[1 4],1); % ref Davidson et.al. 2006
for icell = 1:length(hp)
    
    figure; hold on        
    for iarm = 1:3
        plot(find(armposindex(:,iarm)),sp(armposindex(:,iarm),icell)./ra(armposindex(:,iarm)),'color',newcol(iarm,:),'LineWidth',3)
    end
    title(num2str(hp(icell)))
    axis tight
    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmReplayTriggered\' label '_' thisdir(1:end-4) '_Precision_Raw_HPCell' num2str(hp(icell))])
    toadd = 0;
    figure; hold on        
    for iarm = 1:3
        dat = sp(armposindex(:,iarm),icell)./ra(armposindex(:,iarm));
%         FR=fillmissing(dat,'linear');
%         FR2 = [FR(1)*ones(20,1);FR;FR(end)*ones(20,1)];        
%         FRnew2 = filter2(Filter,FR2(end:-1:1)'); FRnew2 = FRnew2(end:-1:1);
%         FRnew3 = filter2(Filter,FR2');        
%         FRnew4 = mean([FRnew2' FRnew3'],2);
%         smdat=FRnew4(21:end-20);
%         if rem(length(dat),2)==0
%             smdat = nanmean([dat(1:2:end) dat(2:2:end)],2);
%         else
%             smdat = nanmean([dat(1:2:end-1) dat(2:2:end)],2);
%         end
        if rem(length(dat),3)==0
            smdat = nanmean([dat(1:3:end) dat(2:3:end) dat(3:3:end)],2);
        else
            smdat = nanmean([dat(1:3:end-rem(length(dat),3)) dat(2:3:end-rem(length(dat),3)-1) dat(3:3:end)],2);
        end
        
        
        plot([1:length(smdat)]+toadd,smdat,'color',newcol(iarm,:),'LineWidth',3)
        toadd = length(smdat)+toadd;
    end
    axis tight
    title(num2str(hp(icell)))
    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmReplayTriggered\' label '_' thisdir(1:end-4) '_Precision_Binned3_HPCell' num2str(hp(icell))])
end

% %%
%  figure; hold on    
%     
%     for iarm = 1:3
%         plot(find(armposindex(:,iarm)),ra(armposindex(:,iarm)),'color',newcol(iarm,:),'LineWidth',3)
%     end
%     %%
%      figure; hold on    
%     
%     for iarm = 1:3
%         plot(find(armposindex(:,iarm)),sp(armposindex(:,iarm),icell),'color',newcol(iarm,:),'LineWidth',3)
%     end
%         title(num2str(hp(icell)))