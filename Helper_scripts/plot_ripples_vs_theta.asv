% % 'RP_SSDarm2'  %armreplays
% 'RP_Cand_sig_modu_include' %ripples
%  p_SSDreplay %armreplays vs ripples (25% sig, some had very few ripples that
%  arent armevents with these lax criterion)
% Theta modulation


% replaylab = {'AllArmEvents';'Wc3ArmCov3';'WC4ArmCov5mj7';'WC5ArmCov5mj4'};
% ilab = 1; igroup = 1;
% grouplabels = {['AllCells_' replaylab{ilab}];['BestTTDays_' replaylab{ilab}];['BestDays_' replaylab{ilab}]};
% savefolder = ['F:\XY_matdata\Figures\ForPaper\' grouplabels{igroup}];

function plot_ripples_vs_theta(dirs,savefolder,toplot)

close all
cd(dirs.homedir)
d2 = dir('*.mat');
allcells = []; hp_all = [];
for id = 1:size(d2,1)
    thisdir = d2(id).name;             
    [mrv_pval_zval_meandir_kappa,hp_m] = pfc_theta_lock(thisdir,savefolder,toplot);
    load(thisdir, 'RP_Cand_sig_modu_include','other_cells','RP_pSSDreplay')
    if id==1 %(C2)
        examplecell = other_cells==91;
    elseif id==6 %('D3')
        examplecell = other_cells==4;
    elseif id==11 %(I1)
        examplecell = other_cells==83;
    else
        examplecell = false(size(other_cells));
    end
    allcells = [allcells; RP_Cand_sig_modu_include RP_pSSDreplay mrv_pval_zval_meandir_kappa id*ones(size(other_cells)) examplecell];
    hp_all = [hp_all; hp_m];
end

% figure
% a = circ_plot(allcells(allcells(:,1)==1,8),'hist',[],24,true,true,'linewidth',2,'color','r');
% figure
% b = circ_plot(allcells(allcells(:,6)<.05,8),'hist',[],24,true,true,'linewidth',2,'color','r');
% 
% d = mod(rad2deg(circ_mean(allcells(:,8)))+360,360);
% f = mod(rad2deg(circ_median(allcells(:,8)))+360,360);
% e = circ_rtest(allcells(:,8));
% 
% 
% sum(allcells(:,6)<.05)./length(allcells(:,6))




[r,p] = corr(abs(allcells(:,2)),allcells(:,5));
[r2,p2] = corr(log10(abs(allcells(allcells(:,2)~=0,2))),log10(allcells(allcells(:,2)~=0,5)));

if ~isfolder([savefolder '/ThetaBasics/'])
    mkdir([savefolder '/ThetaBasics/'])
end
excells = find(allcells(:,end)==1);

figure;
dat1 = abs(allcells(:,2));
dat2 = allcells(:,5);
loglog(dat1,dat2,'ok')
hold on
for icell = 1:length(excells)    
    loglog(dat1(excells(icell)),dat2(excells(icell)),'o')
end
axis tight
set(gca,'XTick',[.001 .01 .1])
set(gca,'YTick',[.001 .01 .1])
% hold on
% xl = get(gca,'xlim');
% b = polyfit(dat1,dat2,1);
% y2 = polyval(b,xl);
% loglog(xl,y2,'r','LineWidth',3)
set(gcf,'Position',[   2399        -244         215         172])
hold on
%plot bestfit line
set(gca,'FontName','Arial')
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '/ThetaBasics/Ripple_Vs_Theta_LogLog'])

figure; hold on
dat1 = abs(allcells(:,2));
dat2 = allcells(:,5);
plot(dat1,dat2,'ok','MarkerSize',10)
for icell = 1:length(excells)    
    loglog(dat1(excells(icell)),dat2(excells(icell)),'o','MarkerSize',10)
end
axis tight
xl = get(gca,'xlim');
b = polyfit(dat1,dat2,1);
y2 = polyval(b,xl);
plot(xl,y2,'r','LineWidth',3)
%plot bestfit line
xlabel('Ripple modulation')
ylabel('Theta modulation')
set(gcf,'Position',[ 2100        -466         560         420])
yl = get(gca,'ylim');
text(.01,yl(2)*.8,['Linlin: r = ' num2str(r) ' p = ' num2str(p)])
text(.01,yl(2)*.7,['LogLog: r = ' num2str(r2) ' p = ' num2str(p2)])
set(gca,'FontName','Arial')
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '/ThetaBasics/Ripple_Vs_Theta_LinLin'])
% figure; hold on
% plot(log10(abs(allcells(:,2))),log10(allcells(:,5)),'.k')
% set(gcf,'Position',[  2364        -280         215         195])

bboth = sum(allcells(:,6)<.05 & allcells(:,1)==1);
vennX([sum(allcells(:,1)==1)-bboth,bboth,sum(allcells(:,6)<.05)-bboth],.01)
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '/ThetaBasics/RippleModVenn'])

figure; hold on
bb1 = bar([1 2],[sum(allcells(:,1)==1)./size(allcells,1) sum(allcells(:,1)==0)./size(allcells,1); sum(allcells(allcells(:,1)==1,2)>0)./sum(allcells(:,1)==1) sum(allcells(allcells(:,1)==1,2)<0)./sum(allcells(:,1)==1)],'stacked');
th = [sum(allcells(allcells(:,1)==1,6)<.05)./length(allcells(allcells(:,1)==1,6)) sum(allcells(allcells(:,1)==0,6)<.05)./length(allcells(allcells(:,1)==0,6))]*100;
text(0,((sum(allcells(:,1)==1)./size(allcells,1))/2)-.1,[num2str(round(th(1))) '% theta mod'])
text(0,(1-((sum(allcells(:,1)==0)./size(allcells,1))/2))-.1,[num2str(round(th(2))) '% theta mod'])
text(.65,(sum(allcells(:,1)==1)./size(allcells,1))/2,[num2str(sum(allcells(:,1)==1)) ' cells, ' num2str(round(100*sum(allcells(:,1)==1)./size(allcells,1)),2) '%'])
text(.65,1-((sum(allcells(:,1)==0)./size(allcells,1))/2),[num2str(sum(allcells(:,1)==0)) ' cells, ' num2str(round(100*sum(allcells(:,1)==0)./size(allcells,1)),2) '%'])
text(1.65,(sum(allcells(allcells(:,1)==1,2)>0)./sum(allcells(:,1)==1))/2,[num2str(sum(allcells(allcells(:,1)==1,2)>0)) ' Ext, ' num2str(round(100*sum(allcells(allcells(:,1)==1,2)>0)./sum(allcells(:,1)==1)),2) '%'])
text(1.65,1-(sum(allcells(allcells(:,1)==1,2)<0)./sum(allcells(:,1)==1))/2,[num2str(sum(allcells(allcells(:,1)==1,2)<0)) ' Inh, ' num2str(round(100*sum(allcells(allcells(:,1)==1,2)<0)./sum(allcells(:,1)==1)),2) '%'])
legend('Significantly Ripple Modulated','Not Modulated','Location','bestoutside') %,'Positive Modulation','Negative Modulation'
set(gca,'FontName','Arial')
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '/ThetaBasics/RippleMod'])

figure; hold on
bb1 = bar([1 2],[sum(allcells(:,6)<.05)./size(allcells,1) sum(allcells(:,6)>=.05)./size(allcells,1); .5 .5],'stacked');
th = [sum(allcells(allcells(:,6)<.05,1))./length(allcells(allcells(:,6)<.05,1)) sum(allcells(allcells(:,6)>=.05,1))./length(allcells(allcells(:,6)>=.05,1))]*100;
text(0,((sum(allcells(:,6)<.05)./size(allcells,1))/2)-.1,[num2str(round(th(1))) '% ripple mod'])
text(0,(1-((sum(allcells(:,6)>=.05)./size(allcells,1))/2))-.1,[num2str(round(th(2))) '% rpple mod'])
text(.65,(sum(allcells(:,6)<.05)./size(allcells,1))/2,[num2str(sum(allcells(:,6)<.05)) ' cells, ' num2str(round(100*sum(allcells(:,6)<.05)./size(allcells,1)),2) '%'])
text(.65,1-((sum(allcells(:,6)>=.05)./size(allcells,1))/2),[num2str(sum(allcells(:,6)>=.05)) ' cells, ' num2str(round(100*sum(allcells(:,6)>=.05)./size(allcells,1)),2) '%'])
legend('Significantly Theta Modulated','Not Modulated','Location','bestoutside')
set(gca,'FontName','Arial')
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '/ThetaBasics/ThetaMod'])

% [mean(allcells(allcells(:,2)>0 & allcells(:,1)==1,9)) mean(allcells(allcells(:,2)<0 & allcells(:,1)==1,9)) mean(allcells(allcells(:,1)==0,9)) mean(allcells(allcells(:,1)==1,9))]
% [mean(allcells(allcells(:,2)>0 & allcells(:,1)==1,5)) mean(allcells(allcells(:,2)<0 & allcells(:,1)==1,5)) mean(allcells(allcells(:,1)==0,5)) mean(allcells(allcells(:,1)==1,5))]

figure; hold on;
subplot(1,2,1)
aa = circ_plot(hp_m,'hist',[],24,true,true,'linewidth',2,'color','r');
title('HP')
subplot(1,2,2)
c = circ_plot(allcells(:,8),'hist',[],24,true,true,'linewidth',2,'color','r');
title('PFC')
set(gca,'FontName','Arial')
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '/ThetaBasics/Theta_Lock_10cmsec'])


