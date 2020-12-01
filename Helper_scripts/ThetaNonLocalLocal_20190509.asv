cd(dirs.homedir)
d2 = dir('*.mat');
% pcells = []; pcells1 = []; pcells2 = []; celldat = []; shuffdat = []; label = 'RP';
armsig = []; 
numshuff = 1000;
cutoff = .34;
sscells = []; ssScells = [];
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'times_armon_thetaof_headingarm_lap_thetahalf_all','spikedata','other_cells','RP_CandEventTimes','pos','PFCthetaspikes_binned');
    Th = times_armon_thetaof_headingarm_lap_thetahalf_all;
    clear times_armon_thetaof_headingarm_lap_thetahalf_all
    %same as for replay    
%     modind1 = [0 .1];
    bn = squeeze(mean(PFCthetaspikes_binned(:,101:106,:),2));
    
    touse = Th(:,9)>cutoff & ~isnan(Th(:,4));
    Th2 = Th(touse,:);
    bn2 = bn(:,touse);
    clear Th bn


    excludearm = false(3,1);
    for iarm = 1:3        
        IND = Th2(:,3)==iarm;
        uq = unique(Th2(IND,4));
        h = hist(Th2(IND,4),uq);
        numperm = (factorial(sum(h))/(factorial(sum(h)-min(h))*factorial(min(h))));
        if isinf(factorial(sum(h))); numperm = factorial(min(h)); end
        if numperm<numshuff || length(uq)==1
            disp(['Skipping Arm ' num2str(iarm) ' Theta, only ' num2str(numperm) ' perms**'])
            excludearm(iarm,1) = 1;
            continue
        end        
    end
    
    %shuffle for local, keep laps together, shuffle the label of the whole lap           
    L = [Th2(:,3) Th2(:,6)]; % 
    [a,b] = unique(L(:,2));
    laparms = L(b,1);
    LS = NaN(size(L,1),numshuff);
    for ish = 1:numshuff
        laparmsS = laparms(randperm(length(laparms)));        
        for iarm = 1:3            
            LS(ismember(L(:,2),a(laparmsS==iarm)),ish) = iarm;
        end   
        
        while sum(sum(LS(:,ish)~=L(:,1)))==0
            laparmsS = laparms(randperm(length(laparms)));        
            for iarm = 1:3            
                LS(ismember(L(:,2),a(laparmsS==iarm)),ish) = iarm;
            end   
        end
        
        if ish>1
            while any(sum(LS(:,ish)~=LS(:,1:ish-1))==0)
                laparmsS = laparms(randperm(length(laparms)));        
                for iarm = 1:3            
                    LS(ismember(L(:,2),a(laparmsS==iarm)),ish) = iarm;
                end   
            end
        end
    end
    locshuffind = LS;

    %shuffle for local, 
%     L = Th2(:,3); 
%     LS = NaN(size(L,1),numshuff);
%     
%     for ishuff = 1:numshuff
%         las = L(randperm(length(L)));
%         while sum(las~=L)==0
%             las = L(randperm(length(L)));
%         end
%         if ishuff>1
%             while any(sum(las~=LS(:,1:ishuff-1))==0)
%                 las = L(randperm(length(L)));
%             end
%         end
%         LS(:,ishuff) = las;
%     end       
%     locshuffind = LS;
    
    % shuffle for nonlocal, within each arm, shuffle the labels
    L2 = Th2(:,4); 
    LS = NaN(size(L,1),numshuff);
    for iarm = 1:3
        armind = L(:,1)==iarm;
        la = L2(armind);
        if excludearm(iarm)==1
            continue
        end        
        for ishuff = 1:numshuff
            las = la(randperm(length(la)));
            while sum(las~=L(armind,1))==0
                las = la(randperm(length(la)));
            end
            if ishuff>1
                while any(sum(las~=LS(armind,1:ishuff-1))==0)
                    las = la(randperm(length(la)));
                end
            end
            LS(armind,ishuff) = las;
        end        
    end
    nonlocshuffind = LS;
    
    %shuffle nonlocal
%     L = Th2(:,4); 
%     LS = NaN(size(L,1),numshuff);
%     for ishuff = 1:numshuff
%         las = L(randperm(length(L)));
%         while sum(las~=L)==0
%             las = L(randperm(length(L)));
%         end
%         if ishuff>1
%             while any(sum(las~=LS(:,1:ishuff-1))==0)
%                 las = L(randperm(length(L)));
%             end
%         end
%         LS(:,ishuff) = las;
%     end
%     nonlocshuffind = LS;
    
    bnbig = repmat(bn2,[1 1 numshuff]);
    thcell = NaN(size(bn2,1),3,2);  
    thcellshuff = NaN(size(bn2,1),numshuff,3,2);  
    L = Th2(:,3); L2 = Th2(:,4);
    for iarm = 1:3
        if excludearm(iarm)==1
            continue
        end
       nonlocind = Th2(:,4)==iarm;
       locind = Th2(:,3)==iarm;

       thcell(:,iarm,1) = nanmean(bn2(:,L(:,1)==iarm),2);
       thcell(:,iarm,2) = nanmean(bn2(:,L2==iarm),2);

       for ishuff = 1:numshuff
           thcellshuff(:,ishuff,iarm,1) = nanmean(bn2(:,locshuffind(:,ishuff)==iarm),2);
           thcellshuff(:,ishuff,iarm,2) = nanmean(bn2(:,nonlocshuffind(:,ishuff)==iarm),2);
       end           
    end
    
    
    meandat = thcell;
    meandatS = thcellshuff;
    
    cb = nchoosek(1:3,2);
    ss = zeros(size(meandat,1),2);
    ssS = zeros(size(meandat,1),numshuff,2);
    for icomb = 1:3
        ss(:,1) = nansum([ss(:,1) (meandat(:,cb(icomb,1),1)-meandat(:,cb(icomb,2),1)).^2],2);
        ss(:,2) = nansum([ss(:,2) (meandat(:,cb(icomb,1),2)-meandat(:,cb(icomb,2),2)).^2],2);
        for ishuff = 1:numshuff
           ssS(:,ishuff,1) = nansum([ssS(:,ishuff,1) (meandatS(:,ishuff,cb(icomb,1),1)-meandatS(:,ishuff,cb(icomb,2),1)).^2],2);
           ssS(:,ishuff,2) = nansum([ssS(:,ishuff,2) (meandatS(:,ishuff,cb(icomb,1),2)-meandatS(:,ishuff,cb(icomb,2),2)).^2],2);
        end
    end
    sscells = cat(1,sscells,ss);
    ssScells = cat(1,ssScells,ssS);
    load(thisdir,[label '_SSDarm'],[label '_moduarm'],[label '_moduarm2'],[label '_pSSDarm'],[label '_SSDarm2'],[label '_Cand_sig_modu_include'])
    eval(['aa = [' label '_pSSDarm ' label '_SSDarm ' label '_Cand_sig_modu_include ' label '_SSDarm2];'])
    armsig = cat(1,armsig,aa);
    p1 = (sum(nanmean(ssScells(:,:,1))>=nanmean(sscells(:,1)))+1)./(numshuff+1);
    p2 = (sum(nanmean(ssScells(:,:,2))>=nanmean(sscells(:,2)))+1)./(numshuff+1);
    disp([num2str(id) ' ' num2str(p1) ' ' num2str(p2)])    
end

%%
% savelab = 'sigonly';
% ind = armsig(:,1)<.05;
savelab = 'allcells';
ind = true(size(armsig,1),1);
figure; 
histogram(nanmean(ssScells(ind,:,1)),'FaceColor','k'); 
hold on;
histogram(nanmean(ssScells(ind,:,2)),'FaceColor','r'); 

hold on; yl = get(gca,'ylim'); 
plot([nanmean(sscells(ind,1)) nanmean(sscells(ind,1))],yl,'k','LineWidth',3);
plot([nanmean(sscells(ind,2)) nanmean(sscells(ind,2))],yl,'r','LineWidth',3);
helper_savefig(['E:\XY_matdata\Figures\ForPaper\Theta\ThetaNonLocalLocal_pfc_SummaryFigs1_' num2str(cutoff) '_' savelab])

figure; hold on
histogram(nanmean(sscells(ind,1))-nanmean(ssScells(ind,:,1)),'FaceColor','k')
histogram(nanmean(sscells(ind,2))-nanmean(ssScells(ind,:,2)),'FaceColor','r')
helper_savefig(['E:\XY_matdata\Figures\ForPaper\Theta\ThetaNonLocalLocal_pfc_SummaryFigs2_' num2str(cutoff) '_' savelab])
%%
figure; hold on
h1 = histc(realdat,-6:.5:6); 
h2 = histc(nanmean(shuff),-.5:.05:.5); 
% b = bar([h1./sum(h1) h2./sum(h2)],'LineWidth',2,'BarWidth',1);
% b(2).FaceColor = [0 0 0];
% b(1).FaceColor = [1 1 1];
% bar(-6:.1:6,h1./sum(h1),'FaceColor','k','FaceAlpha',.5,'EdgeAlpha',1)
bar(-.5:.05:.5,h2./sum(h2),'FaceColor','w','FaceAlpha',.5,'LineWidth',2)
yl = get(gca,'ylim');
plot([nanmean(realdat) nanmean(realdat)],yl,'k','LineWidth',3)
xlabel('Firing Rate Difference Between Theta Sweeps')
ylabel('Count')
set(gca,'FontName','Arial','FontSize',18)

axes('Position',[.16 .65 .2 .2])
box on
hold on
pe = pie([sum(~isnan(pcells2)) sum(pcells2<.05)]);
pe(1).FaceColor = 'w';
pe(3).FaceColor = 'k';
p= 1-binocdf(sum(pcells2<.05)-1,sum(~isnan(pcells2)),.05);
title(['Cells Significantly Modulated' newline ' by Direction of Theta Sweep'])
xlim([-1.7 1.7])
ylim([-1.5 1.5])
set(gca,'xtick',[],'ytick',[])

set(gcf,'Position',[  680   388   816   590])

set(gca,'FontName','Arial')

set(gcf,'renderer','Painters')
helper_savefig(['E:\XY_matdata\Figures\ForPaper\Theta\NonLocalThetaDifference_pfc_SummaryFigs3_' num2str(cutoff) '_' savelab])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Figure3\NonLocalThetaDifference_pfc_SummaryFigs3_' num2str(cutoff) '_' savelab])






