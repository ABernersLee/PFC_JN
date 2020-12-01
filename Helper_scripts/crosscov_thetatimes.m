function crosscov_thetatimes(dirs,toplot,toplotcells,tosave,speedcutoff,globe,igroup,savefolder)
cd(dirs.homedir)

d2 = dir('*.mat');


if ~isfolder([savefolder '\Theta\ThetaXcorr\'])
    mkdir([savefolder '\Theta\ThetaXcorr\'])
end
[FilterA] =fspecial('gaussian',[5 1],5/3);
thetacycdiff = NaN(size(d2,1),1);
crossesall = []; crossesraw = [];
load(['F:\XY_matdata\AllDays\matfiles\theta_sequence_zero.mat'],'theta_sequence_zero')
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'vel','hp_cells','other_cells','hpinterneurons','spikedata','pos','theta_globalzero','HP_Theta','other_cells_touse')
    other_cells = other_cells(other_cells_touse(:,igroup));
    
    theta =  HP_Theta; %HPTheta_spikephase;
    
    if globe == 1
        savelab = 'GlobalZero';
        theta(:,4) = mod(theta(:,4)-theta_globalzero,360);
    else
        savelab = 'SequenceZero_extrapolated';
        theta(:,4) = mod((theta(:,4)-theta_globalzero)-theta_sequence_zero,360);
    end

    hpcells = hp_cells(~ismember(hp_cells,hpinterneurons));
    clear hp_cells hpinterneurons

    % only periods of time where velocity is over the speed cutoff for at least
    % a second
    abv = vel>speedcutoff & diff([pos(:,1);pos(end,1)])<5; %  %going above speed cutoff
    dabv = [diff([0;abv])];
    st = pos(dabv==1,1);
    nd = pos(dabv==-1,1);
    st = st(1:length(nd));
    tms = [st nd];
    tms(diff(tms,[],2)<1,:) = [];

    s = spikedata;
%     tms = tms';    
%     [~,i] = histc(spikes(:,1),tms(:));
%     spikes(rem(i,2)==0,:) = []; %taking out spikes that aren't in these periods

    
    %for armsig_mx_SequenceZero3 or 4 way (downsampling)
%     [~,u] = histc(s(:,1),theta(:,1));
%     s(u==0,:) = [];
%     u(u==0) = [];
%     s(:,4) = theta(u,4);
    
    %for armsig_mx_SequenceZero2 way (extrapolating)
    s2 = [s NaN(size(s,1),2)]; theta1 = theta;
    ind = find(theta1(:,4)>358 & theta1(:,4)<360);
    ind(diff(ind)==1) = [];
    anglechange = theta1(ind,1);
    as = [anglechange(1:end-1) anglechange(2:end)];
    ind2 = (as(:,2)-as(:,1))<.1 | (as(:,2)-as(:,1))>.2;    
    as(ind2,:) = [];  
    spikeind = find(s(:,1)>as(1,1) & s(:,1)<as(end,1));    
    ss = s(spikeind,:);
    for it = 1:size(as,1)
        spikeind2 = ss(:,1)>=as(it,1) & ss(:,1)<as(it,2);
        s3 = ss(spikeind2,:);
        cycleind = spikeind(spikeind2);
        a = (s3(:,1)-as(it,1))/(as(it,2)-as(it,1));
        s2(cycleind,end-1:end) = [ones(size(s3,1),1)*it mod(a*360,360)]; 
        clear s3 a
    end
    s(:,4) = s2(:,5);

    earlyt = theta(:,4)<180; latet = theta(:,4)>180;    
    thetacycdiff(id,1) = 100*((sum(latet)-sum(earlyt))./sum(earlyt));
    
    %this is the armsig_mx_SequenceZero3 and 4 way (downsampling)
%     if sum(latet)>sum(earlyt)
%         ffind = find(latet);
%         rind = randperm(length(ffind));        
%         theta(ffind(rind(1:sum(latet)-sum(earlyt))),4) = NaN;
%     elseif sum(earlyt)>sum(latet)
%         ffind = find(earlyt);
%         rind = randperm(length(ffind));        
%         theta(ffind(rind(1:sum(earlyt)-sum(latet))),4) = NaN;
%     end
    
    earlyt = theta(:,4)<180; latet = theta(:,4)>180;
    if sum(earlyt)~=sum(latet)
        disp('did not work')
    end
    
    figure; histogram(range(tms,2),20)
    helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\' num2str(thisdir(1:end-4)) '_timerange'])
    
    bin = .01;
    wind = .25;
    % crosses = NaN(round((.125/bin))*2+1,length(hpcells),length(other_cells),3);
    crosses = NaN(round((wind/bin))*2+1,length(other_cells),3); %until     9/18/2020
    crosses1 = NaN(round((wind/bin))*2+1,size(tms,1),length(other_cells),3); %made 9/18/2020
    tic
    
    
    for itms = 1:size(tms,1) %added 9/18/2020 instead of zeroing out bins that I take out, do each period seperately
        tmstime = s(:,1)>= tms(itms,1) & s(:,1)<= tms(itms,2);
        T = range(tms(itms,:));

%         %way one %pairs of hp-pfc neurons        
%         h1 = NaN(length(tms(itms,1):bin:tms(itms,2)),length(hpcells),3);
%         for ihcell = 1:length(hpcells)
%             h1(:,ihcell,1) = histc(s(s(:,2)==hpcells(ihcell) & tmstime,1),tms(itms,1):bin:tms(itms,2));
%             h1(:,ihcell,2) = histc(s(s(:,2)==hpcells(ihcell) & s(:,4)<180 & tmstime,1),tms(itms,1):bin:tms(itms,2));
%             h1(:,ihcell,3) = histc(s(s(:,2)==hpcells(ihcell) & s(:,4)>180 & tmstime,1),tms(itms,1):bin:tms(itms,2));
%         end
% 
%         %way one %pairs of hp-pfc neurons
%         for ipfc = 1:length(other_cells)
%             p1 = histc(s(s(:,2)==other_cells(ipfc) & tmstime,1),tms(itms,1):bin:tms(itms,2));  
%             if size(p1,1)==1; p1 = p1'; end
%             for ih = 1:3
%                 p2 = [p1 h1(:,:,ih)];
%                 [cc,lags] = xcov(p2,round((wind/bin)),'coeff');
%                 cct = cc(:,2:size(p2,2));            
%                 crosses1(:,itms,ipfc,ih) = nanmean(cct,2) * T; 
%             end
%         end        
% % %way two with pairs (not the right nidex yet but would be faster)
% %         p1 = NaN(length(tms(itms,1):bin:tms(itms,2)),length(other_cells));
% %         for ipfc = 1:length(other_cells)
% %             p1(:,ipfc) = histc(s(s(:,2)==other_cells(ipfc) & tmstime,1),tms(itms,1):bin:tms(itms,2));  
% %         end
% %         for ih = 1:3
% %             p2 = [p1 h1(:,:,ih)];            
% %             [cc,lags] = xcov(p2,round((wind/bin)),'coeff');
% %             Nlab = repmat(1:size(p2,2),[1 size(p2,2)]);
% %             Nlab2 = [1:size(p2,2)]'*ones(1,size(p2,2));
% %             Nlab2 = Nlab2(:);
% %             cct = cc(:,Nlab>1 & Nlab2'<length(other_cells)); %this isn't             always 1
% %             cct2 = reshape(cct,[size(cct,1) length(other_cells) length(hpcells)]);
% %             crosses1(:,itms,ipfc,ih) = nanmean(cct,2) * T; 
% %         end

        %bulk hp activity (looks the same as pairs)
        p1 = NaN(length(tms(itms,1):bin:tms(itms,2)),length(other_cells));
        for ipfc = 1:length(other_cells)
            p1(:,ipfc) = histc(s(s(:,2)==other_cells(ipfc) & tmstime,1),tms(itms,1):bin:tms(itms,2));  
        end 
        h1 = NaN(length(tms(itms,1):bin:tms(itms,2)),3);
        h1(:,1) = histc(s(ismember(s(:,2),hpcells) & tmstime,1),tms(itms,1):bin:tms(itms,2));
        h1(:,2) = histc(s(ismember(s(:,2),hpcells) & s(:,4)<180 & tmstime,1),tms(itms,1):bin:tms(itms,2));
        h1(:,3) = histc(s(ismember(s(:,2),hpcells) & s(:,4)>180 & tmstime,1),tms(itms,1):bin:tms(itms,2));
        for ih = 1:3
            p2 = [h1(:,ih) p1];
            [cc,lags] = xcov(p2,round((wind/bin)),'coeff');            
            crosses1(:,itms,:,ih) = cc(end:-1:1,2:size(p2,2)) * sqrt(bin*T); %opposite order as used to be (hp cell, then pfc cell)
        end
        

        if rem(itms,50)==0
            ttoc = toc;
            disp([num2str(itms) ' of ' num2str(size(tms,1)) ' done (' num2str(round(100*itms/size(tms,1))) '%) in ' num2str(ttoc/60) ' min'])
        end
    end
    

    c = squeeze(nanmean(crosses1,2));        
    for ih = 1:3
        h1 = filtfilt(FilterA,1,[repmat(c(1,:,ih),[5 1 1]); c(:,:,ih); repmat(c(end,:,ih),[5 1 1])]);                       
        crosses(:,:,ih) = h1(6:end-5,:);
    end  

    if toplotcells
          for ipfc = 1:length(other_cells)
            figure; hold on;
            lagsp = lags*bin;                

%             fwd = [crosses(:,ipfc,1)-semcross(:,1)]'; rev = [crosses(:,ipfc,1)+semcross(:,1)]';
%             p2 = patch([lagsp lagsp(end:-1:1)],[fwd rev(end:-1:1)],'black');
%             p2.FaceAlpha=.1; p2.EdgeAlpha=0;
%             fwd = [crosses(:,ipfc,2)-semcross(:,2)]'; rev = [crosses(:,ipfc,2)+semcross(:,2)]';
%             p2 = patch([lagsp lagsp(end:-1:1)],[fwd rev(end:-1:1)],'blue');
%             p2.FaceAlpha=.1; p2.EdgeAlpha=0;
%             fwd = [crosses(:,ipfc,3)-semcross(:,3)]'; rev = [crosses(:,ipfc,3)+semcross(:,3)]';
%             p2 = patch([lagsp lagsp(end:-1:1)],[fwd rev(end:-1:1)],'red');
%             p2.FaceAlpha=.1; p2.EdgeAlpha=0;

            a = plot(lagsp,crosses(:,ipfc,1),'k','LineWidth',3);
            aa = plot(lagsp,crosses(:,ipfc,2),'b','LineWidth',3);
            aaa = plot(lagsp,crosses(:,ipfc,3),'r','LineWidth',3);

            yl = get(gca,'ylim');
            plot([0 0],yl,'k','LineWidth',3)
            title(['PFC cell ' num2str(other_cells(ipfc))])
            xlabel('Seconds')
            xlim([-wind wind])
            ylabel('Cross-Covariance')
            set(gca,'FontSize',18,'ylim',yl)
            legend([a aa aaa],{'All HP Spikes';'First Half Theta';'Second Half Theta'})
            set(gcf,'Position',[ 680         121        1098         857])
            helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\' savelab '_' num2str(thisdir(1:end-4)) '_PFCcell' num2str(other_cells(ipfc)) '_newtms'])
          end
    end
     
    if toplot
    
        figure; hold on;
%         crosses2 = cat(1,crosses(:,:,2),crosses(:,:,3));
%         crosses3 = zscore(crosses2);
%         crossesplot = cat(3,crosses3(1:size(crosses,1),:),crosses3(size(crosses,1)+1:size(crosses,1)*2,:));
%         crossesplot = zscore(crosses(:,:,2:3));
        crossesplot = crosses(:,:,2:3);
        semcross = squeeze(nanstd(crossesplot,[],2))./repmat(sqrt(squeeze(sum(sum(~isnan(crossesplot))>0)))',[size(crossesplot,1) 1]);
        lagsp = lags*bin;
        
        aa = plot(lagsp,nanmean(crossesplot(:,:,1),2),'b','LineWidth',3);
        fwd = [nanmean(crossesplot(:,:,1),2)-semcross(:,1)]'; rev = [nanmean(crossesplot(:,:,1),2)+semcross(:,1)]';
        p2 = patch([lagsp lagsp(end:-1:1)],[fwd rev(end:-1:1)],'blue');
        p2.FaceAlpha=.1; p2.EdgeAlpha=0;
        
        aaa = plot(lagsp,nanmean(crossesplot(:,:,2),2),'r','LineWidth',3);
        fwd = [nanmean(crossesplot(:,:,2),2)-semcross(:,2)]'; rev = [nanmean(crossesplot(:,:,2),2)+semcross(:,2)]';
        p2 = patch([lagsp lagsp(end:-1:1)],[fwd rev(end:-1:1)],'red');
        p2.FaceAlpha=.1; p2.EdgeAlpha=0;
        
        yl = get(gca,'ylim');
        plot([0 0],yl,'k','LineWidth',3)    
        xlabel('Seconds')
        xlim([-wind wind])
        ylabel('Cross-Covariance')
        set(gca,'FontSize',18,'ylim',yl)
        title(savelab)
        legend([aa aaa],{'First Half Theta';'Second Half Theta'},'Location','northwest')
        set(gcf,'Position',[ 680         121        1098         857])
        helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\avgPFCcells_' num2str(thisdir(1:end-4)) '_' savelab '_firstsecondonly_new'])
    end
    
    if tosave
        
%         crosses_raw = c;
%         [mx1a,mx1] = max(crosses_raw(round(wind/bin)+1:end,:,2));
%         [mx1a(2,:),mx1(2,:)] = max(crosses_raw(round(wind/bin)+1:end,:,3));
        
        [mx1a,mx1] = max(crosses(round(wind/bin)+1:end,:,2));
        [mx1a(2,:),mx1(2,:)] = max(crosses(round(wind/bin)+1:end,:,3));

        [~,mx2] = max(mx1a);
        mx2(sum(isnan(mx1))==2) = NaN;

        label = 'RP';
        load(thisdir, [label '_moduarm'],[label '_pSSDarm'],[label '_SSDarm'],[label '_Cand_sig_modu_include'])
        eval(['armsig = [' label '_pSSDarm ' label '_SSDarm ' label '_Cand_sig_modu_include];'])
        armsig_mx = armsig;
        armsig_mx(:,6) = mx2;
        armsig_mx(:,7:8) = mx1a';
        
        %new ones as of 9/4/2020 after reviewer comments (2-3)
        %new one with tms on 9/18/2020 after review comments (4)
        if globe == 1
            armsig_mx_GlobalZero2 = armsig_mx;
            save(thisdir,'armsig_mx_GlobalZero2','-append')        
        else
            armsig_mx_SequenceZero2 = armsig_mx;
            save(thisdir,'armsig_mx_SequenceZero2','-append')       
        end
    end
    
    
    disp([num2str(itms) ' of ' num2str(size(tms,1)) ' done (' num2str(round(100*itms/size(tms,1))) '%) in ' num2str(ttoc/60) ' min'])
    crossesall = cat(2,crossesall,crosses);
    crossesraw = cat(2,crossesraw,c);
    disp([num2str(id)])
end

ssd = [];
 for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'RP_SSDarm2','armsig_mx_SequenceZero4')
    ssd = cat(1,ssd,[armsig_mx_SequenceZero4 RP_SSDarm2]);
 end
 
    lagsall2 = NaN(size(crossesall,2),2);
    lags2 = lags(round(wind/bin)+1:find(lags*bin<.125,1,'last'));
    for ipfc = 1:size(crossesall,2)
        lagsall2(ipfc,1) = lags2(max(crossesall(round(wind/bin)+1:find(lags*bin<.125,1,'last'),ipfc,1))==(crossesall(round(wind/bin)+1:find(lags*bin<.125,1,'last'),ipfc,1)))*bin;
        lagsall2(ipfc,2) = lags2(max(crossesall(round(wind/bin)+1:find(lags*bin<.125,1,'last'),ipfc,3))==(crossesall(round(wind/bin)+1:find(lags*bin<.125,1,'last'),ipfc,3)))*bin;
    end
    
%      lagsall2 = NaN(size(crossesall,2),2);
%     lags2 = lags(round(wind/bin)+1:end);
%     for ipfc = 1:size(crossesall,2)
%         lagsall2(ipfc,1) = lags2(max(crossesall(round(wind/bin)+1:end,ipfc,1))==(crossesall(round(wind/bin)+1:end,ipfc,1)))*bin;
%         lagsall2(ipfc,2) = lags2(max(crossesall(round(wind/bin)+1:end,ipfc,3))==(crossesall(round(wind/bin)+1:end,ipfc,3)))*bin;
%     end
    
if toplot
%     [FilterA] =fspecial('gaussian',[5 1],5/3);
%     crossZ = zscore(crossesall);    
%     
%     figure; hold on;
%     for ih = 1:3
%         h1 = filter(FilterA,1,nanmean([repmat(crossZ(1,:,ih),[5 1 1]); crossZ(:,:,ih); repmat(crossZ(end,:,ih),[5 1 1])],2)*bin*2);
%         plot(lags*bin,h1(6:end-5));
%     end
%     yl = get(gca,'ylim');
%     plot([0 0], yl, 'k--')
%     plot([.06 .06], yl, 'k--')
%     legend({'Both';'FirstHalf';'SecondHalf'})
    
    lagsp = lags*bin;
    
     for ih = 1:3
        figure; hold on;

        crossZ = zscore(crossesall(:,ssd(:,3)==1,ih),[],1);   
%         h1 = filter(FilterA,1,nanmean([repmat(crossZ(1,:),[5 1 ]); crossZ; repmat(crossZ(end,:),[5 1 ])],2));
%         sem1 = nanstd(crossZ,[],2)./sqrt(sum(sum(~isnan(crossZ))>0));
%         sem = filter(FilterA,1,[repmat(sem1(1),[5 1]); sem1; repmat(sem1(end),[5 1])]);
%         
        h1 = nanmean(crossZ,2);
        sem = nanstd(crossZ,[],2)./sqrt(sum(sum(~isnan(crossZ))>0));
        

        fwd = [h1-sem]'; rev = [h1+sem]';
        p2 = patch([lagsp lagsp(end:-1:1)],[fwd rev(end:-1:1)],'black');
        p2.FaceAlpha=.1; p2.EdgeAlpha=0;
        plot(lagsp,h1,'k','LineWidth',3);
        yl = get(gca,'ylim');
        plot([0 0], yl, 'k--')
        plot([.06 .06], yl, 'k--')
         xlabel('Seconds')
        xlim([-wind wind])
        ylabel('Cross-Covariance (zscored)')
        set(gca,'FontSize',18,'ylim',yl)
        title(savelab)
        helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\SigRippleMod_vgPFCcellsZscore_Delay_' savelab '_' num2str(ih)])

        figure; hold on;
        crossZ = crossesall(:,ssd(:,3)==1,ih);  
        h1 = nanmean(crossZ,2);
        sem = nanstd(crossZ,[],2)./sqrt(sum(sum(~isnan(crossZ))>0));

        fwd = [h1-sem]'; rev = [h1+sem]';
        p2 = patch([lagsp lagsp(end:-1:1)],[fwd rev(end:-1:1)],'black');
        p2.FaceAlpha=.1; p2.EdgeAlpha=0;
        plot(lagsp,h1,'k','LineWidth',3);
        yl = get(gca,'ylim');
        plot([0 0], yl, 'k--')
        plot([.06 .06], yl, 'k--')
         xlabel('Seconds')
        xlim([-wind wind])
        ylabel('Cross-Covariance')
        set(gca,'FontSize',18,'ylim',yl)
        title(savelab)
        helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\SigRippleMod_PFCcells_Delay_' savelab '_' num2str(ih)])
     end
    
    for ih = 1:3
        figure; hold on;

        crossZ = zscore(crossesall(:,:,ih));   
%         h1 = filter(FilterA,1,nanmean([repmat(crossZ(1,:),[5 1 ]); crossZ; repmat(crossZ(end,:),[5 1 ])],2));
%         sem1 = nanstd(crossZ,[],2)./sqrt(sum(sum(~isnan(crossZ))>0));
%         sem = filter(FilterA,1,[repmat(sem1(1),[5 1]); sem1; repmat(sem1(end),[5 1])]);
%         
        h1 = nanmean(crossZ,2);
        sem = nanstd(crossZ,[],2)./sqrt(sum(sum(~isnan(crossZ))>0));
        

        fwd = [h1-sem]'; rev = [h1+sem]';
        p2 = patch([lagsp lagsp(end:-1:1)],[fwd rev(end:-1:1)],'black');
        p2.FaceAlpha=.1; p2.EdgeAlpha=0;
        plot(lagsp,h1,'k','LineWidth',3);
        yl = get(gca,'ylim');
        plot([0 0], yl, 'k--')
        plot([.06 .06], yl, 'k--')
         xlabel('Seconds')
        xlim([-wind wind])
        ylabel('Cross-Covariance (zscored)')
        set(gca,'FontSize',18,'ylim',yl)
        title(savelab)
        helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\avgPFCcellsZscore_Delay_' savelab '_' num2str(ih)])

        figure; hold on;
        crossZ = crossesall(:,:,ih);  
        h1 = nanmean(crossZ,2);
        sem = nanstd(crossZ,[],2)./sqrt(sum(sum(~isnan(crossZ))>0));

        fwd = [h1-sem]'; rev = [h1+sem]';
        p2 = patch([lagsp lagsp(end:-1:1)],[fwd rev(end:-1:1)],'black');
        p2.FaceAlpha=.1; p2.EdgeAlpha=0;
        plot(lagsp,h1,'k','LineWidth',3);
        yl = get(gca,'ylim');
        plot([0 0], yl, 'k--')
        plot([.06 .06], yl, 'k--')
         xlabel('Seconds')
        xlim([-wind wind])
        ylabel('Cross-Covariance')
        set(gca,'FontSize',18,'ylim',yl)
        title(savelab)
        helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\avgPFCcells_Delay_' savelab '_' num2str(ih)])
    end
    
    figure; hold on;    
    crosses2 = cat(1,crossesall(:,:,1),crossesall(:,:,2),crossesall(:,:,3));
    crosses3 = zscore(crosses2);
    crossesplot = cat(3,crosses3(1:size(crossesall,1),:),crosses3(size(crossesall,1)+1:size(crossesall,1)*2,:),crosses3(size(crossesall,1)*2+1:end,:));
    plot(lags*bin,nanmean(crossesplot(:,:,1),2),'k','LineWidth',3)
    plot(lags*bin,nanmean(crossesplot(:,:,2),2),'b','LineWidth',3)
    plot(lags*bin,nanmean(crossesplot(:,:,3),2),'r','LineWidth',3)
    yl = get(gca,'ylim');
    plot([0 0],yl,'k','LineWidth',3)    
    xlabel('Seconds')
    xlim([-wind wind])
    ylabel('Cross-Covariance (Z-scored)')
    set(gca,'FontSize',18,'ylim',yl)
    title(savelab)
    legend({'All HP Spikes';'First Half Theta';'Second Half Theta'},'Location','northwest')
    set(gcf,'Position',[ 680         121        1098         857])
    helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\avgPFCcells_Zscore_' savelab])
    
    figure; hold on;    
    crosses2 = cat(1,crossesall(:,:,1),crossesall(:,:,2),crossesall(:,:,3));
    crosses3 = crosses2;
    crossesplot = cat(3,crosses3(1:size(crossesall,1),:),crosses3(size(crossesall,1)+1:size(crossesall,1)*2,:),crosses3(size(crossesall,1)*2+1:end,:));
    plot(lags*bin,nanmean(crossesplot(:,:,1),2),'k','LineWidth',3)
    plot(lags*bin,nanmean(crossesplot(:,:,2),2),'b','LineWidth',3)
    plot(lags*bin,nanmean(crossesplot(:,:,3),2),'r','LineWidth',3)
    yl = get(gca,'ylim');
    plot([0 0],yl,'k','LineWidth',3)    
    xlabel('Seconds')
    xlim([-wind wind])
    ylabel('Cross-Covariance')
    set(gca,'FontSize',18,'ylim',yl)
    title(savelab)
    legend({'All HP Spikes';'First Half Theta';'Second Half Theta'},'Location','northwest')
    set(gcf,'Position',[ 680         121        1098         857])
    helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\avgPFCcells_' savelab])
    
        figure; hold on;    
    crosses2 = cat(1,crossesall(:,ssd(:,3)==1,1),crossesall(:,ssd(:,3)==1,2),crossesall(:,ssd(:,3)==1,3));
    crosses3 = zscore(crosses2);
    crossesplot = cat(3,crosses3(1:size(crossesall,1),:),crosses3(size(crossesall,1)+1:size(crossesall,1)*2,:),crosses3(size(crossesall,1)*2+1:end,:));
    plot(lags*bin,nanmean(crossesplot(:,:,1),2),'k','LineWidth',3)
    plot(lags*bin,nanmean(crossesplot(:,:,2),2),'b','LineWidth',3)
    plot(lags*bin,nanmean(crossesplot(:,:,3),2),'r','LineWidth',3)
    yl = get(gca,'ylim');
    plot([0 0],yl,'k','LineWidth',3)    
    xlabel('Seconds')
    xlim([-wind wind])
    ylabel('Cross-Covariance (Z-scored)')
    set(gca,'FontSize',18,'ylim',yl)
    title(savelab)
    legend({'All HP Spikes';'First Half Theta';'Second Half Theta'},'Location','northwest')
    set(gcf,'Position',[ 680         121        1098         857])
    helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\SigRippleMod_PFCcells_Zscore_' savelab])
    
    figure; hold on;    
    crosses2 = cat(1,crossesall(:,ssd(:,3)==1,1),crossesall(:,ssd(:,3)==1,2),crossesall(:,ssd(:,3)==1,3));
    crosses3 = crosses2;
    crossesplot = cat(3,crosses3(1:size(crossesall,1),:),crosses3(size(crossesall,1)+1:size(crossesall,1)*2,:),crosses3(size(crossesall,1)*2+1:end,:));
    plot(lags*bin,nanmean(crossesplot(:,:,1),2),'k','LineWidth',3)
    plot(lags*bin,nanmean(crossesplot(:,:,2),2),'b','LineWidth',3)
    plot(lags*bin,nanmean(crossesplot(:,:,3),2),'r','LineWidth',3)
    yl = get(gca,'ylim');
    plot([0 0],yl,'k','LineWidth',3)    
    xlabel('Seconds')
    xlim([-wind wind])
    ylabel('Cross-Covariance')
    set(gca,'FontSize',18,'ylim',yl)
    title(savelab)
    legend({'All HP Spikes';'First Half Theta';'Second Half Theta'},'Location','northwest')
    set(gcf,'Position',[ 680         121        1098         857])
    helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\SigRippleMod_PFCcells_' savelab])
    
%     figure; hold on;
%     crosses2 = cat(1,crossesall(:,:,2),crossesall(:,:,3));
%     crosses3 = zscore(crosses2);
%     crossesplot = cat(3,crosses3(1:size(crossesall,1),:),crosses3(size(crossesall,1)+1:size(crossesall,1)*2,:));
%     
%     plot(lags*bin,nanmean(crossesplot(:,:,1),2),'b','LineWidth',3)
%     plot(lags*bin,nanmean(crossesplot(:,:,2),2),'r','LineWidth',3)
%     yl = get(gca,'ylim');
%     plot([0 0],yl,'k','LineWidth',3)    
%     xlabel('Seconds')
%     xlim([-wind wind])
%     ylabel('Cross-Covariance')
%     set(gca,'FontSize',18,'ylim',yl)
%     title(savelab)
%     legend({'First Half Theta';'Second Half Theta'},'Location','northwest')
%     set(gcf,'Position',[ 680         121        1098         857])
%     helper_saveandclosefig([savefolder '\Theta\ThetaXcorr\avgPFCcells_' savelab '_firstsecondonly_new'])
end