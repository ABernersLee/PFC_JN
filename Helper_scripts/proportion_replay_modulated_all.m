function proportion_replay_modulated_all(label,savelab,igroup,savefolder)

cd('F:\XY_matdata\AllDays')

d2 = dir('*.mat');

%%
SMI = []; PARM = []; 
dday = [];
for id =  1:size(d2,1)
        thisdir = d2(id).name;
        load(thisdir,[label '_Cand_sig_modu_include'],[label '_pSSDarm'],'other_cells_touse')
        eval(['smi = ' label '_Cand_sig_modu_include;']);
        eval(['parm = ' label '_pSSDarm;']);        
        
        SMI = cat(1,SMI,smi(other_cells_touse(:,igroup),:));
        PARM = cat(1,PARM,parm(other_cells_touse(:,igroup),:));
        dday = cat(1,dday,id*ones(size(other_cells_touse,1),1));
        clear smi parm
end
%%        
ss = SMI;
% ss(SMI(:,3)~=1,1:2) = NaN; 
% PARM(SMI(:,3)~=1,1) = NaN;
%%


figure; hold on
b = [sum(ss(:,1)==0)./sum(~isnan(ss(:,1))) sum(ss(:,1)==1)./sum(~isnan(ss(:,1))) ; ...    
    sum(ss(ss(:,1)==1,2)<0)./sum(ss(:,1)==1) sum(ss(ss(:,1)==1,2)>0)./sum(ss(:,1)==1) ; ...
    sum(PARM(:,1)>=.05)./sum(~isnan(PARM(:,1))) sum(PARM(:,1)<.05)./sum(~isnan(PARM(:,1))) ; ...
     sum(ss(PARM(:,1)<.05,2)<0)./sum(PARM(:,1)<.05) sum(ss(PARM(:,1)<.05,2)>0)./sum(PARM(:,1)<.05)];
bb1 = bar([1 4],b([1 3],end:-1:1),'stacked');
bb1(1).BarWidth = .2; bb1(2).BarWidth = .2;
bb1(2).FaceColor = [.7 .7 .7]; %[.5 .5 .5]; 
bb1(1).FaceColor = 'none'; %'k'
bb2 = bar([2 5],b([2 4],:),'stacked');
bb2(1).BarWidth = .2; bb2(2).BarWidth = .2;
bb2(1).FaceColor = 'b'; bb2(2).FaceColor = 'r';
xlim([0 6])
set(gca,'ytick',0:.2:1,'yticklabel',0:20:100)
ylabel('Percent')
set(gca,'xtick',[1.5 4.5],'xticklabel',{'Ripple Triggered';'Arm-Replay Modulated'})


text(.85,.9,num2str(round(b(1,1)*100)),'Color','k','FontSize',18);
text(.8,.2,num2str(round(b(1,2)*100)),'Color','k','FontSize',18);

text(1.85,.2,num2str(round(b(2,1)*100)),'Color','w','FontSize',18);
text(1.85,.9,num2str(round(b(2,2)*100)),'Color','w','FontSize',18);

plot([3.5 4.5], [.05 .05],'r--')
legend('Not Modulated','Significantly Modulated','Negative Modulation','Positive Modulation','Chance','Location','bestoutside')
text(3.7,.9,[num2str(round(b(3,1)*100)) '%'],'Color','k','FontSize',18);
text(3.75,.85,num2str(sum(PARM(:,1)>=.05)),'Color','k','FontSize',18);
text(3.7,.15,[num2str(round(b(3,2)*100)) '%'],'Color','k','FontSize',18);
text(3.85,.1,num2str(sum(PARM(:,1)<.05)),'Color','k','FontSize',18);


text(4.85,.2,num2str(round(b(4,1)*100)),'Color','w','FontSize',18);
text(4.85,.9,num2str(round(b(4,2)*100)),'Color','w','FontSize',18);

set(gca,'FontSize',20)
set(gcf,'Position',[2072         225         995         609])
if ~isfolder([savefolder '\ReplayTriggeredStats\'])
    mkdir([savefolder '\ReplayTriggeredStats\'])
end
helper_saveandclosefig([savefolder '\ReplayTriggeredStats\' label '_BarCharts_' savelab])
%%
figure; hold on


subplot(3,2,1)
a = sum(ss(:,1)==0);
b = sum(ss(:,1)==1);
p1 = (1-binocdf(b-1,a,.05));
pp = pie([a b],{['Not Modulated (' num2str(sum(ss(:,1)==0)) ')'];['Modulated (' num2str(sum(ss(:,1)==1)) '), p = ' num2str(round(p1,2,'significant'))]});
pp(1).FaceColor = 'k';
pp(3).FaceColor = [.5 .5 .5];
if p1<.05
    pp(4).Color = 'r';
    pp(3).FaceColor = 'r';
end

title(['Proportion of Cells Modulated by ' label])


subplot(3,2,2)
a = sum(ss(ss(:,1)==1,2)<0);
b = sum(ss(ss(:,1)==1,2)>0);
p1 = (1-binocdf(b-1,a,.05));
pp = pie([a b],{['Negativly Modulated (' num2str(sum(ss(ss(:,1)==1,2)<0)) ')'];['Positively Modulated (' num2str(sum(ss(ss(:,1)==1,2)>0)) '), p = ' num2str(round(p1,2,'significant'))]});
pp(1).FaceColor = 'k';
pp(3).FaceColor = [.5 .5 .5];
if p1<.05
    pp(4).Color = 'r';
    pp(3).FaceColor = 'r';
end

title(['Direction of Modulation of Sig Modulated Cells'])

subplot(3,2,3)
a = sum(PARM(:,1)>=.05);
b = sum(PARM(:,1)<.05);
p1 = (1-binocdf(b-1,a,.05));
pp = pie([a b],{['Not Modulated (' num2str(sum(PARM(:,1)>=.05)) ')'];['Modulated (' num2str(sum(PARM(:,1)<.05)) '), p = ' num2str(round(p1,2,'significant'))]});
title(['Proportion of Cells Modulated by Specific Arms in Replays'])
pp(1).FaceColor = 'k';
pp(3).FaceColor = [.5 .5 .5];
if p1<.05
    pp(4).Color = 'r';
    pp(3).FaceColor = 'r';
end


subplot(3,2,4)
a = sum(ss(PARM(:,1)<.05,2)<0);
b = sum(ss(PARM(:,1)<.05,2)>0);
p1 = (1-binocdf(b-1,a,.05));
pp = pie([a b],{['Negativly Modulated (' num2str(sum(ss(PARM(:,1)<.05,2)<0)) ')'];['Positively Modulated (' num2str(sum(ss(PARM(:,1)<.05,2)>0)) '), p = ' num2str(round(p1,2,'significant'))]});
title(['Direction of Modulation of Arm-Replay Modulated Cells'])
pp(1).FaceColor = 'k';
if length(pp)>2
    pp(3).FaceColor = [.5 .5 .5];
    if p1<.05
        pp(4).Color = 'r';
        pp(3).FaceColor = 'r';
    end
end

subplot(3,2,5)
a = sum(PARM(ss(:,1)==1,1)>=0.05);
b = sum(PARM(ss(:,1)==1,1)<.05);
p1 = (1-binocdf(b-1,a,.05));
pp = pie([a b],{['Not Modulated (' num2str(sum(PARM(ss(:,1)==1,1)>=0.05)) ')'];['Modulated (' num2str(sum(PARM(ss(:,1)==1,1)<.05)) '), p = ' num2str(round(p1,2,'significant'))]});
title(['Proportion of ' label ' Mod Cells also Arm-Replay Modulated'])
pp(1).FaceColor = 'k';
pp(3).FaceColor = [.5 .5 .5];
if p1<.05
    pp(4).Color = 'r';
    pp(3).FaceColor = 'r';
end

set(gcf,'Position',[1993         -53        1126        1032])


helper_saveandclosefig([savefolder '\ReplayTriggeredStats\' label '_PieCharts_' savelab])

    %%
    
            
figure; hold on
ss = SMI;
ss(SMI(:,3)~=1,1:2) = NaN; 
PARM(SMI(:,3)~=1,1) = NaN;

subplot(3,2,1)
a = sum(ss(:,1)==0);
b = sum(ss(:,1)==1);
p1 = (1-binocdf(b-1,a,.05));
pp = pie([a b],{['Not Modulated (' num2str(sum(ss(:,1)==0)) ')'];['Modulated (' num2str(sum(ss(:,1)==1)) '), p = ' num2str(round(p1,2,'significant'))]});
pp(1).FaceColor = 'k';
pp(3).FaceColor = [.5 .5 .5];
if p1<.05
    pp(4).Color = 'r';
    pp(3).FaceColor = 'r';
end

title(['Proportion of Cells Modulated by ' label])


subplot(3,2,2)
a = sum(ss(ss(:,1)==1,2)<0);
b = sum(ss(ss(:,1)==1,2)>0);
p1 = (1-binocdf(b-1,a,.05));
pp = pie([a b],{['Negativly Modulated (' num2str(sum(ss(ss(:,1)==1,2)<0)) ')'];['Positively Modulated (' num2str(sum(ss(ss(:,1)==1,2)>0)) '), p = ' num2str(round(p1,2,'significant'))]});
pp(1).FaceColor = 'k';
pp(3).FaceColor = [.5 .5 .5];
if p1<.05
    pp(4).Color = 'r';
    pp(3).FaceColor = 'r';
end

title(['Direction of Modulation of Sig Modulated Cells'])

subplot(3,2,3)
a = sum(PARM(:,1)>=.05);
b = sum(PARM(:,1)<.05);
p1 = (1-binocdf(b-1,a,.05));
pp = pie([a b],{['Not Modulated (' num2str(sum(PARM(:,1)>=.05)) ')'];['Modulated (' num2str(sum(PARM(:,1)<.05)) '), p = ' num2str(round(p1,2,'significant'))]});
title(['Proportion of Cells Modulated by Specific Arms in Replays'])
pp(1).FaceColor = 'k';
pp(3).FaceColor = [.5 .5 .5];
if p1<.05
    pp(4).Color = 'r';
    pp(3).FaceColor = 'r';
end


subplot(3,2,4)
a = sum(ss(PARM(:,1)<.05,2)<0);
b = sum(ss(PARM(:,1)<.05,2)>0);
p1 = (1-binocdf(b-1,a,.05));
pp = pie([a b],{['Negativly Modulated (' num2str(sum(ss(PARM(:,1)<.05,2)<0)) ')'];['Positively Modulated (' num2str(sum(ss(PARM(:,1)<.05,2)>0)) '), p = ' num2str(round(p1,2,'significant'))]});
title(['Direction of Modulation of Arm-Replay Modulated Cells'])
pp(1).FaceColor = 'k';
if length(pp)>2
    pp(3).FaceColor = [.5 .5 .5];
    if p1<.05
        pp(4).Color = 'r';
        pp(3).FaceColor = 'r';
    end
end

subplot(3,2,5)
a = sum(PARM(ss(:,1)==1,1)>=0.05);
b = sum(PARM(ss(:,1)==1,1)<.05);
p1 = (1-binocdf(b-1,a,.05));
pp = pie([a b],{['Not Modulated (' num2str(sum(PARM(ss(:,1)==1,1)>=0.05)) ')'];['Modulated (' num2str(sum(PARM(ss(:,1)==1,1)<.05)) '), p = ' num2str(round(p1,2,'significant'))]});
title(['Proportion of ' label ' Mod Cells also Arm-Replay Modulated'])
pp(1).FaceColor = 'k';
pp(3).FaceColor = [.5 .5 .5];
if p1<.05
    pp(4).Color = 'r';
    pp(3).FaceColor = 'r';
end

set(gcf,'Position',[1993         -53        1126        1032])

        
helper_saveandclosefig([savefolder '\ReplayTriggeredStats\' label '_PieCharts_ExludinglowFR' savelab])


figure; hold on
b = [sum(ss(:,1)==0)./sum(~isnan(ss(:,1))) sum(ss(:,1)==1)./sum(~isnan(ss(:,1))) ; ...    
    sum(ss(ss(:,1)==1,2)<0)./sum(ss(:,1)==1) sum(ss(ss(:,1)==1,2)>0)./sum(ss(:,1)==1) ; ...
    sum(PARM(:,1)>=.05)./sum(~isnan(PARM(:,1))) sum(PARM(:,1)<.05)./sum(~isnan(PARM(:,1))) ; ...
     sum(ss(PARM(:,1)<.05,2)<0)./sum(PARM(:,1)<.05) sum(ss(PARM(:,1)<.05,2)>0)./sum(PARM(:,1)<.05)];
bb1 = bar([1 4],b([1 3],:),'stacked');
bb1(1).BarWidth = .2; bb1(2).BarWidth = .2;
bb1(1).FaceColor = [.5 .5 .5]; bb1(2).FaceColor = 'k';
bb2 = bar([2 5],b([2 4],:),'stacked');
bb2(1).BarWidth = .2; bb2(2).BarWidth = .2;
bb2(1).FaceColor = 'b'; bb2(2).FaceColor = 'r';
xlim([0 6])
set(gca,'ytick',0:.2:1,'yticklabel',0:20:100)
ylabel('Percent')
set(gca,'xtick',[1.5 4.5],'xticklabel',{'Ripple Triggered';'Arm-Replay Modulated'})
legend('Not Modulated','Significantly Modulated','Negative Modulation','Positive Modulation','Location','bestoutside')

text(.85,.2,num2str(round(b(1,1)*100)),'Color','w','FontSize',18);
text(.8,.9,num2str(round(b(1,2)*100)),'Color','w','FontSize',18);

text(1.85,.2,num2str(round(b(2,1)*100)),'Color','w','FontSize',18);
text(1.85,.9,num2str(round(b(2,2)*100)),'Color','w','FontSize',18);

text(3.85,.2,num2str(round(b(3,1)*100)),'Color','w','FontSize',18);
text(3.85,.9,num2str(round(b(3,2)*100)),'Color','w','FontSize',18);

text(4.85,.2,num2str(round(b(4,1)*100)),'Color','w','FontSize',18);
text(4.85,.9,num2str(round(b(4,2)*100)),'Color','w','FontSize',18);

set(gca,'FontSize',20)
set(gcf,'Position',[2072         225         995         609])
helper_saveandclosefig([savefolder '\ReplayTriggeredStats\' label '_BarCharts_ExcludinglowFR' savelab])

%% from the ones that already are ripple modulated


ss = SMI;

figure; hold on
b = [sum(ss(:,1)==0)./sum(~isnan(ss(:,1))) sum(ss(:,1)==1)./sum(~isnan(ss(:,1))) ; ...    
    sum(ss(ss(:,1)==1,2)<0)./sum(ss(:,1)==1) sum(ss(ss(:,1)==1,2)>0)./sum(ss(:,1)==1) ; ...
    sum(PARM(:,1)>=.05 & ss(:,1)==1)./sum(ss(:,1)==1 & ~isnan(PARM(:,1))) sum(PARM(:,1)<.05 & ss(:,1)==1)./sum(ss(:,1)==1 & ~isnan(PARM(:,1))) ; ...
     sum(ss(PARM(:,1)<.05 & ss(:,1)==1,2)<0)./sum(PARM(:,1)<.05 & ss(:,1)==1) sum(ss(PARM(:,1)<.05 & ss(:,1)==1,2)>0)./sum(PARM(:,1)<.05 & ss(:,1)==1)];
bb1 = bar([1 4],b([1 3],end:-1:1),'stacked');
bb1(1).BarWidth = .2; bb1(2).BarWidth = .2;
bb1(1).FaceColor = 'none'; % [.5 .5 .5];
bb1(2).FaceColor = 'none';
bb2 = bar([2 5],b([2 4],:),'stacked');
bb2(1).BarWidth = .2; bb2(2).BarWidth = .2;
bb2(1).FaceColor = 'b'; bb2(2).FaceColor = 'r';
xlim([0 6])
set(gca,'ytick',0:.2:1,'yticklabel',0:20:100)
ylabel('Percent')
set(gca,'xtick',[1.5 4.5],'xticklabel',{'Ripple Triggered';'Arm-Replay Modulated'})
legend('Not Modulated','Significantly Modulated','Negative Modulation','Positive Modulation','Location','bestoutside')


text(.85,.9,num2str(round(b(1,1)*100)),'Color','k','FontSize',18);
text(.8,.2,num2str(round(b(1,2)*100)),'Color','k','FontSize',18);

text(1.85,.2,num2str(round(b(2,1)*100)),'Color','w','FontSize',18);
text(1.85,.9,num2str(round(b(2,2)*100)),'Color','w','FontSize',18);

plot([3.5 4.5], [.05 .05],'r--')
text(3.7,.9,[num2str(round(b(3,1)*100)) '%'],'Color','k','FontSize',18);
text(3.75,.85,num2str( sum(PARM(:,1)>=.05 & ss(:,1)==1)),'Color','k','FontSize',18);
text(3.7,.15,[num2str(round(b(3,2)*100)) '%'],'Color','k','FontSize',18);
text(3.85,.1,num2str( sum(PARM(:,1)<.05 & ss(:,1)==1)),'Color','k','FontSize',18);


text(4.85,.2,num2str(round(b(4,1)*100)),'Color','w','FontSize',18);
text(4.85,.9,num2str(round(b(4,2)*100)),'Color','w','FontSize',18);

set(gca,'FontSize',20)
set(gcf,'Position',[2072         225         995         609])
helper_saveandclosefig([savefolder '\ReplayTriggeredStats\' label '_RippleModBarCharts_' savelab])

%%
ss(SMI(:,3)~=1,1:2) = NaN; 
PARM(SMI(:,3)~=1,1) = NaN;

figure; hold on
b = [sum(ss(:,1)==0)./sum(~isnan(ss(:,1))) sum(ss(:,1)==1)./sum(~isnan(ss(:,1))) ; ...    
    sum(ss(ss(:,1)==1,2)<0)./sum(ss(:,1)==1) sum(ss(ss(:,1)==1,2)>0)./sum(ss(:,1)==1) ; ...
    sum(PARM(:,1)>=.05 & ss(:,1)==1)./sum(ss(:,1)==1) sum(PARM(:,1)<.05 & ss(:,1)==1)./sum(ss(:,1)==1) ; ...
     sum(ss(PARM(:,1)<.05 & ss(:,1)==1,2)<0)./sum(PARM(:,1)<.05 & ss(:,1)==1) sum(ss(PARM(:,1)<.05 & ss(:,1)==1,2)>0)./sum(PARM(:,1)<.05 & ss(:,1)==1)];

bb1 = bar([1 4],b([1 3],:),'stacked');
bb1(1).BarWidth = .2; bb1(2).BarWidth = .2;
bb1(1).FaceColor = [.5 .5 .5]; bb1(2).FaceColor = 'k';
bb2 = bar([2 5],b([2 4],:),'stacked');
bb2(1).BarWidth = .2; bb2(2).BarWidth = .2;
bb2(1).FaceColor = 'b'; bb2(2).FaceColor = 'r';
xlim([0 6])
set(gca,'ytick',0:.2:1,'yticklabel',0:20:100)
ylabel('Percent')
set(gca,'xtick',[1.5 4.5],'xticklabel',{'Ripple Triggered';'Arm-Replay Modulated'})
legend('Not Modulated','Significantly Modulated','Negative Modulation','Positive Modulation','Location','bestoutside')

text(.85,.2,num2str(round(b(1,1)*100)),'Color','w','FontSize',18);
text(.8,.9,num2str(round(b(1,2)*100)),'Color','w','FontSize',18);

text(1.85,.2,num2str(round(b(2,1)*100)),'Color','w','FontSize',18);
text(1.85,.9,num2str(round(b(2,2)*100)),'Color','w','FontSize',18);

text(3.85,.2,num2str(round(b(3,1)*100)),'Color','w','FontSize',18);
text(3.85,.9,num2str(round(b(3,2)*100)),'Color','w','FontSize',18);

text(4.85,.2,num2str(round(b(4,1)*100)),'Color','w','FontSize',18);
text(4.85,.9,num2str(round(b(4,2)*100)),'Color','w','FontSize',18);

set(gca,'FontSize',20)
set(gcf,'Position',[2072         225         995         609])
helper_saveandclosefig([savefolder '\ReplayTriggeredStats\' label '_RippleModBarCharts_ExcludinglowFR' savelab])

% %%
% mmlab = {'All Cells';'Neg Mod';'Pos Mod';'Sig Neg Mod';'Sig Pos Mod'};
% rplab = {'Ripple Cand Events';'SD Cand Events'};
% for irpsd = 1:2
%     if irpsd == 1
%         dat = allcells1;
%     elseif irpsd == 2
%         dat = allcells2;
%     end
%     for mm = 4:5    
%         if mm == 1
%             ind = smi(:,3)==1;
%         elseif mm == 2 %neg mod
%             ind = smi(:,3)==1 & smi(:,2)<0;
%         elseif mm == 3 %pos mod
%             ind = smi(:,3)==1 & smi(:,2)>0;
%         elseif mm == 4 % sig neg mod
%             ind = smi(:,3)==1 & smi(:,2)<0 & smi(:,1)==1;
%         elseif mm==5 %sig pos mod
%             ind = smi(:,3)==1 & smi(:,2)>0 & smi(:,1)==1;
%         end
% 
%         r1 = corr([1:size(dat,2)]',dat(ind,:)');
%         nanmean(r1);
%         figure; hold on
%         sem = nanstd(zscore(dat(ind,:)'),[],2)./sqrt(size(dat(ind,:),1));
%         meandat = nanmean(zscore(dat(ind,:)'),2);
%         rev = meandat-sem;
%         plot(meandat,'k','LineWidth',2)
%         patch([1:size(dat,2) size(dat,2):-1:1],[meandat+sem;rev(end:-1:1)],'black','FaceAlpha',.3,'EdgeAlpha',0)
%         [r2,p2] = corr([1:size(dat,2)]',meandat);
%         tt = title([mmlab{mm} ' ' rplab{irpsd} ' r = ' num2str(r2) ' p = ' num2str(p2)]);
%         if p2<.05
%             tt.Color = 'r';
%         end
%     end
% end