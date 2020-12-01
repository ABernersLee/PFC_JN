function Figure1_PlaceFields(dirs,savefolder)
cd(dirs.homedir)
d2 = dir('*.mat');
hpf = [];
pfcf = [];
binsize = 20;
SigmaField = 2;
Two_D_Filter=fspecial('gaussian',[10 10],SigmaField);
% Two_D_Filter=fspecial('gaussian',[3 3],SigmaField);
velcutoff = 5;
fieldfigs = false;
tocalc = false;
for id = 1:size(d2,1)
    clear OpenFR
    load(d2(id).name,'InFR','OutFR','spikedata','vel','hp_cells','other_cells','hpinterneurons','pos','OpenFR','armposindex')
    hpf = cat(1,hpf,mean(cat(3,InFR(hp_cells(~ismember(hp_cells,hpinterneurons)),:),OutFR(hp_cells(~ismember(hp_cells,hpinterneurons)),:)),3));
    pfcf = cat(1,pfcf,mean(cat(3,InFR(other_cells,:),OutFR(other_cells,:)),3));        
    if tocalc
        pos2 = pos(:,2:3);
        pos2(:,1) = ceil(pos2(:,1)/binsize);
        pos2(:,2) = ceil(pos2(:,2)/binsize);
        pos2 = pos2-min(pos2)+1;
        tm = diff(pos(:,1)); tm = [tm;tm(end)];

        pos3 = sub2ind(max(pos2),pos2(:,1),pos2(:,2));
        pos3(vel<velcutoff) = NaN;
        ind = 1:max(pos2(:,1))*max(pos2(:,2));
        occ = NaN(size(ind,2),1);
        f = zeros(size(ind,2),max(spikedata(:,2)));
        for j = 1:length(ind)
            occ(ind(j)) = sum(tm(pos3==ind(j)));
            [h,k] = histc(spikedata(pos3(spikedata(:,3)) == ind(j),2),1:max(spikedata(:,2)));
            if ~isempty(k)
                f(ind(j),k) = h(k);            
            end
            if rem(j,1000)==0
                disp(['j = ' num2str(j)])
            end
        end
        occ2 = reshape(occ,[max(pos2)]);        
        f2 = reshape(f,[max(pos2) size(f,2)]);

        OpenFR = f2./repmat(occ2,[1 1 size(f2,3)]);
        save(d2(id).name,'OpenFR','-append')
    else
        load(d2(id).name,'OpenFR')
    end
        
    if fieldfigs        
%         if exist('OpenFR','var'); tocalc = false; else; tocalc = true; end
                
        figure; hold on
        dat = OpenFR(:,:,1);
        imagesc(dat); 
%         colormap autumn
        colormap gray
        cm = get(gca,'colormap');
%         cm = cm(end:-1:1,:);
        close gcf

        cm2 = [cm(end-1:-1:1,:); [1 1 1]];
    
        for icell = 1:max(spikedata(:,2))
            if ismember(icell,other_cells); cellid = 'PFC'; elseif ismember(icell,hpinterneurons); cellid = 'HP INT'; elseif ismember(icell,hp_cells); cellid = 'HP'; end
            figure; hold on
            dat = OpenFR(:,:,icell);
            datz = reshape(nanzscore(dat(:)),[size(dat,1) size(dat,2)]);
            datnan = isnan(datz);
            datz(datnan) = 0;
            dat = filter2(Two_D_Filter,datz);
            smx = max(max(dat));
            dat(datnan) = smx+3*range(range(dat))/size(cm2,1);
            imagesc(dat); 
            set(gca,'colormap',cm2)
            axis xy
            cmi = get(gca,'clim');
%             set(gca,'clim',[0 cmi(2)])
            colorbar('Ticks',[min(min(dat)) max(max(dat))],'TickLabels',[round(min(min(dat)),2,'significant') round(max(max(dat)),2,'significant')])
            title([d2(id).name(1:end-4) ', Cell ' num2str(icell) ' ' cellid])
            helper_saveandclosefig([savefolder '\Figure1\bwZ_' cellid '_' d2(id).name(1:end-4) '_Cell' num2str(icell)]) 
        end
        
        figure; hold on        
        datnan = isnan(dat);
        imagesc(datnan); 
        set(gca,'colormap',cm)
        axis xy
        if ~isfolder([savefolder '\Figure1\IndividualCells\'])
            mkdir([savefolder '\Figure1\IndividualCells\'])
        end
        helper_saveandclosefig([savefolder '\Figure1\IndividualCells\bw_' d2(id).name(1:end-4)])
        
    end
        
end
    
%%
arm11 = pfcf(:,armposindex(:,1));
arm1 = mean(cat(3,arm11(:,1:2:end-1),arm11(:,2:2:end)),3); clear arm11
arm2 = pfcf(:,armposindex(:,2));
arm3 = pfcf(:,armposindex(:,3));
arm2 = arm2(:,1:end-1);
arm3 = arm3(:,1:end-1);
arms = cat(3,arm1,arm2,arm3);
rp = NaN(size(arms,1),3);
for icell = 1:size(arms,1)
    jnk = corr(squeeze(arms(icell,:,:)));
    rp(icell,:) = [jnk(1,2:3) jnk(2,3)];
end
r_pfc = mean(rp,2);
%%

arm11 = hpf(:,armposindex(:,1));
arm1 = mean(cat(3,arm11(:,1:2:end-1),arm11(:,2:2:end)),3); clear arm11
arm2 = hpf(:,armposindex(:,2));
arm3 = hpf(:,armposindex(:,3));
arm2 = arm2(:,1:end-1);
arm3 = arm3(:,1:end-1);
arms = cat(3,arm1,arm2,arm3);
rp = NaN(size(arms,1),3);
for icell = 1:size(arms,1)
    jnk = corr(squeeze(arms(icell,:,:)));
    rp(icell,:) = [jnk(1,2:3) jnk(2,3)];
end
r_hp = mean(rp,2);
%%
figure; hold on
h1 = histcounts(r_hp,-.4:.1:1)'; 
[h2,c] = histcounts(r_pfc,-.4:.1:1); h2 = h2';
% b = bar([h1./sum(h1) h2./sum(h2)],'LineWidth',2,'BarWidth',1);
plot(c(1:end-1)+diff(c)/2,h1./sum(h1),'.-','color',[.5 .5 .5],'LineWidth',2,'MarkerSize',20)
plot(c(1:end-1)+diff(c)/2,h2./sum(h2),'.-','color',[0 0 0],'LineWidth',2,'MarkerSize',20)
% b(2).FaceColor = [0 0 0];
% b(1).FaceColor = [1 1 1];
% set(gca,'xtick',1:2:length(c),'xticklabel',c(1:2:end))
% yl = get(gca,'ylim');
% plot([find(c==0) find(c==0)],yl,'r--','LineWidth',2)
legend('HP','PFC')
ylabel('Probability (Count/Total)')
xlabel('Radial Symmetry (Correlation across arms)')
set(gca,'FontSize',18)
axis tight
set(gcf,'renderer','Painters')
if ~isfolder([savefolder '\Figure1\'])
    mkdir([savefolder '\Figure1\'])
end
[p,h,stats] = ranksum(r_pfc',r_hp');
% [h,p,stats] = ttest2(r_pfc',r_hp');
helper_saveandclosefig([savefolder '\Figure1\radial_prob_nobar'])
% The theta statistic [U/(n1*n2)] is provided as an additional effect size reflecting the distance between the two underlying frequency distributions from which the samples are drawn (Newcombe, 2006a). It is equivalent to the area under the receiver operating characteristic (ROC) curve.
figure; hold on
text(.1,.9,['N = ' num2str(length(r_pfc)) ', ' num2str(length(r_hp))])
text(.1,.8,['p = ' num2str(p)])
text(.1,.7,['U = ' num2str(stats.ranksum)])
text(.1,.6,['theta = ' num2str(stats.ranksum/(length(r_pfc)*length(r_hp)))])
text(.1,.5,['zval = ' num2str(stats.zval)])

helper_saveandclosefig([savefolder '\Figure1\radial_prob_nobar_stats'])

%%
close all
[m2,m] = max(hpf,[],2);
[~,mm] = sort(m);
dat = hpf(mm,:)./repmat(m2(mm),[1 size(hpf,2)]);

figure; imagesc(-dat)
axis xy
colormap gray
helper_saveandclosefig([savefolder '\Figure1\imagesc_HP'])

figure; hold on
iii = 0;
for ii = 3:6:size(dat,1)
    iii = iii+1;
    dat2 = dat(ii,:)+iii-1;
    patch([1:size(dat,2) size(dat,2):-1:1]',[(iii-1)*ones(size(dat2)) dat2(end:-1:1)]','black','EdgeAlpha',0)
end
ylim([0 158])
xlim([1 size(dat,2)])
set(gca,'ytick',[10.5:10:157.5],'yticklabel',[10:10:158])
set(gcf,'Position',[ 2648         -55         517         987])
xlabel('Linearized Position on Track (cm)')
ylabel('HP Neuron')
set(gca,'FontSize',18)
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Figure1\patch_HP'])

[m2,m] = max(pfcf,[],2);
[~,mm] = sort(m);
dat = pfcf(mm,:)./repmat(m2(mm),[1 size(pfcf,2)]);
figure; imagesc(-dat)
axis xy
colormap gray
helper_saveandclosefig([savefolder '\Figure1\imagesc_PFC'])


figure; hold on
iii = 0;
for ii = 1:size(dat,1)
    iii = iii+1;
    dat2 = dat(ii,:)+iii-1;
    patch([1:size(dat,2) size(dat,2):-1:1]',[(iii-1)*ones(size(dat2)) dat2(end:-1:1)]','black')
end

ylim([0 158])
xlim([1 size(dat,2)])
set(gca,'ytick',[10.5:10:157.5],'yticklabel',[10:10:158])
set(gcf,'Position',[ 2648         -55         517         987])
set(gca,'FontSize',18)
xlabel('Linearized Position on Track (cm)')
ylabel('PFC Neuron')
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Figure1\patch_PFC'])

%%
sel_hp = max(hpf,[],2)./mean(hpf,2);
sel_pfc = max(pfcf,[],2)./mean(pfcf,2);
[p,h,stats] = ranksum(sel_pfc',sel_hp');
figure; hold on
text(.1,.9,['N = ' num2str(length(r_pfc)) ', ' num2str(length(r_hp))])
text(.1,.8,['p = ' num2str(p)])
text(.1,.7,['U = ' num2str(stats.ranksum)])
text(.1,.6,['theta = ' num2str(stats.ranksum/(length(r_pfc)*length(r_hp)))])
text(.1,.5,['zval = ' num2str(stats.zval)])
helper_saveandclosefig([savefolder '\Figure1\selectivity_prob2_nobar_stats'])
%%
figure; hold on
histogram(sel_hp,1:37,'FaceColor','k','LineWidth',1); 
hold on; histogram(sel_pfc,1:37,'FaceColor','r','LineWidth',1); 
legend('hp','pfc')
set(gca,'FontSize',18)
ylabel('Count (Neurons)')
xlabel('Spacial Selectivity (Max/Mean)')
set(gca,'FontSize',18)
axis tight
helper_saveandclosefig([savefolder '\Figure1\selectivity_count'])
%%
figure; hold on
% h = histogram(sel_hp,1:37,'Normalization','probability','FaceColor','k','LineWidth',1); 
% hold on; histogram(sel_pfc,1:37,'Normalization','probability','FaceColor','r','LineWidth',1); 
h1 = histc(sel_hp,1:37); 
h2 = histc(sel_pfc,1:37); 
bar(1:37,h1./sum(h1),'FaceColor','k','FaceAlpha',.5,'EdgeAlpha',1)
bar(1:37,h2./sum(h2),'FaceColor','r','FaceAlpha',.5,'EdgeAlpha',1)
legend('hp','pfc')
ylabel('Probability (Count/Total)')
xlabel('Spacial Selectivity (Max/Mean)')
set(gca,'FontSize',18)
axis tight
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Figure1\selectivity_prob'])
%%

figure; hold on

% h = histogram(sel_hp,1:37,'Normalization','probability','LineWidth',1,'DisplayStyle','stairs'); 
% hold on; histogram(sel_pfc,1:37,'Normalization','probability','LineWidth',1,'DisplayStyle','stairs');
h1 = histcounts(sel_hp,1:15)'; 
h2 = histcounts(sel_pfc,1:15)'; 
h1 = [h1; sum(sel_hp>15)];
h2 = [h2; sum(sel_pfc>15)];
% b = bar([h1./sum(h1) h2./sum(h2)],'LineWidth',2,'BarWidth',1);
% b(2).FaceColor = [0 0 0];
% b(1).FaceColor = [1 1 1];
plot(1:length(h2),h1./sum(h1),'.-','color',[.5 .5 .5],'LineWidth',2,'MarkerSize',20)
plot(1:length(h2),h2./sum(h2),'.-','color',[0 0 0],'LineWidth',2,'MarkerSize',20)

set(gca,'xtick',1:2:15,'xticklabel',[1:2:15 15])
% bar(1:37,h1./1sum(h1),'FaceColor','k','FaceAlpha',.5,'EdgeAlpha',1)
% bar(1:37,h2./sum(h2),'FaceColor','r','FaceAlpha',.5,'EdgeAlpha',1)
legend('hp','pfc')
ylabel('Probability (Count/Total)')
xlabel('Spacial Selectivity (Max/Mean)')
set(gca,'FontSize',18)
axis tight
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Figure1\selectivity_prob2_nobar'])