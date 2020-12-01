%older, use Figure3_Internal1d.m

EstBin = .04; cellcutoff = 5; cmcutoff = 5; 
dat = []; datS = []; 
SSD = []; SSDp = [];
toplot = true;
numshuff = 50;
for id = 1:size(d2,1)    
%     if id == 1; toplot = true; else; toplot = false; end
    [SSDarm,alldat,alldatS] = internal_fields(d2(id).name,EstBin,cellcutoff,cmcutoff,toplot,numshuff);
    dat = cat(1,dat,alldat);
    datS = cat(1,datS,alldatS);
    load(d2(id).name,'RP_moduarm','RP_pSSDarm')
    SSD = cat(1,SSD,RP_moduarm(:,RP_pSSDarm<.05)');
%     SSD = cat(1,SSD,RP_moduarm');
%     SSDp = cat(1,SSDp,RP_pSSDarm);
    disp(id)
end


for ilab = 1:3
    
real = squeeze((dat(:,1,:,2,ilab,1)-dat(:,2,:,2,ilab,1)).^2)+...
    squeeze((dat(:,1,:,2,ilab,1)-dat(:,3,:,2,ilab,1)).^2)+...
    squeeze((dat(:,2,:,2,ilab,1)-dat(:,3,:,2,ilab,1)).^2);
shuff = (squeeze(datS(:,1,:,2,ilab,1,:)-datS(:,2,:,2,ilab,1,:)).^2)+...
    squeeze((datS(:,1,:,2,ilab,1,:)-datS(:,3,:,2,ilab,1,:)).^2)+...
    squeeze((datS(:,2,:,2,ilab,1,:)-datS(:,3,:,2,ilab,1,:)).^2);
sigdat = (sum(shuff>=real,3)+1)./(numshuff+1);
disp(sum(sigdat<.05)./size(sigdat,1)*100)

% percent of sig modulated (low velocity, high velocity)
% ilab 1 = 18.5185   25.9259
% 
% ilab 2 = 40.7407   37.0370
% 
% il1b 3 = 44.4444   29.6296

% percent all cells
% 22.1519   22.1519
% 
%    18.9873   14.5570
% 
%    19.6203   18.3544


% dat = NaN(length(SSDarm),3,2,2,3,2); 
%cells,arms,velocity (low, high), locality (local, nonlocal), measurment (mean, median, cmcutoff), raw or smoothed
% raw data across cells

% label = {'Mean','Median',num2str(cmcutoff)};
label = {'Mean_allcells','Median_allcells',[num2str(cmcutoff) '_allcells']};
% label = {'Mean_sigcells','Median_sigcells',[num2str(cmcutoff) '_sigcells']};


%nonlocal high veloccity
figure; hold on
subplot(1,2,1); hold on
dat2 = dat(:,:,2,2,ilab,1);
plot(SSD(:),dat2(:),'ok','MarkerSize',10,'LineWidth',3)
xl = get(gca,'xlim');
xx = SSD(:); yy = dat2(:);
lm = polyfit(xx(~isnan(yy) & ~isnan(xx)),yy(~isnan(yy) & ~isnan(xx)),1);
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',2)
[r,p] = corr(SSD(:),dat2(:),'rows','complete');
xlabel('Change in FR to Replays of Different Arms')
ylabel('FR to Non-Local HP Representation while running')
title(['Non-Local: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))])
set(gca,'FontSize',18)

%local high veloccity
subplot(1,2,2); hold on
dat2 = dat(:,:,2,1,ilab,1);
dat2(isinf(dat2)) = NaN;
plot(SSD(:),dat2(:),'ok','MarkerSize',10,'LineWidth',3)
xl = get(gca,'xlim');
xx = SSD(:); yy = dat2(:);
lm = polyfit(xx(~isnan(yy) & ~isnan(xx)),yy(~isnan(yy) & ~isnan(xx)),1);
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',2)
[r,p] = corr(SSD(:),dat2(:),'rows','complete');
xlabel('Change in FR to Replays of Different Arms')
ylabel('FR to Local HP Representation while running')
title(['Local: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))])
set(gca,'FontSize',18)
suptitle(['Across ' num2str(size(dat2,1)) ' PFC cells'])
set(gcf,'Position',[2126         -40        1351         725])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\AcrossDays_NonLocalVLocal_Raw_' label{ilab}])

%nonlocal high veloccity
% figure; hold on
% dat2 = dat(:,:,2,2,2,1);
% r1 = NaN(size(SSD,1),1); p = r;
% for icell = 1:size(SSD,1)
%     [r1(icell),p(icell)] = corr(SSD(icell,:)',dat2(icell,:)');
% end
% 
% 
% %local high veloccity
% dat2 = dat(:,:,2,1,2,1);
% r2 = NaN(size(SSD,1),1); p = r;
% for icell = 1:size(SSD,1)
%     [r2(icell),p(icell)] = corr(SSD(icell,:)',dat2(icell,:)');
% end
% 
% histogram(r2,5,'FaceColor','k')
% histogram(r1,5,'FaceColor','r')
% p = ranksum(r1,r2,'tail','right');
% title(['Non-Local vs Local: p = ' num2str(p)])
% suptitle('Across PFC cells')

% mean subtract for each cell

%nonlocal high veloccity
figure; hold on
subplot(1,2,1); hold on
dat2 = dat(:,:,2,2,ilab,1);
dat3 = (dat2-median(dat2,2)); % median subtracted
SSD2 = (SSD-median(SSD,2));
% dat3 = zscore(dat2,[],2); %zscored
% SSD2 = zscore(SSD,[],2);
plot(SSD2(:),dat3(:),'ok','MarkerSize',10,'LineWidth',3)
xl = get(gca,'xlim');
xx = SSD2(:); yy = dat3(:);
lm = polyfit(xx(~isnan(yy) & ~isnan(xx)),yy(~isnan(yy) & ~isnan(xx)),1);
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',2)
[r,p] = corr(SSD2(:),dat3(:),'rows','complete','Type','Spearman');
title(['r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))])
xlabel('Change in FR to Replays of Different Arms')
ylabel('FR to Non-Local HP Representation while running')
set(gca,'FontSize',18)
%local high veloccity
subplot(1,2,2); hold on

dat2 = dat(:,:,2,1,ilab,1);
dat2(isinf(dat2)) = NaN;
dat3 = (dat2-median(dat2,2));
SSD2 = (SSD-median(SSD,2));
% dat3 = nanzscore(dat2,[],2); % zscored
% SSD2 = zscore(SSD,[],2);
plot(SSD2(:),dat3(:),'ok','MarkerSize',10,'LineWidth',3)
xl = get(gca,'xlim');
xx = SSD2(:); yy = dat3(:);
lm = polyfit(xx(~isnan(yy) & ~isnan(xx)),yy(~isnan(yy) & ~isnan(xx)),1);
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',2)
[r,p] = corr(SSD2(:),dat3(:),'rows','complete','Type','Spearman');
xlabel('Change in FR to Replays of Different Arms')
ylabel('FR to Local HP Representation while running')
title(['r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))])
set(gca,'FontSize',18)
suptitle('Across PFC cells, Median subtracted from each cell')
set(gcf,'Position',[2126         -40        1351         725])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\AcrossDays_NonLocalVLocal_MedianSubtracted_' label{ilab}])
% zscored

%nonlocal high veloccity
figure; hold on
subplot(1,2,1); hold on
dat2 = dat(:,:,2,2,ilab,1);
% dat3 = (dat2-median(dat2,2)); % median subtracted
% SSD2 = (SSD-median(SSD,2));
dat3 = zscore(dat2,[],2); %zscored
SSD2 = zscore(SSD,[],2);
plot(SSD2(:),dat3(:),'ok','MarkerSize',10,'LineWidth',3)
xl = get(gca,'xlim');
xx = SSD2(:); yy = dat3(:);
lm = polyfit(xx(~isnan(yy) & ~isnan(xx)),yy(~isnan(yy) & ~isnan(xx)),1);
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',2)
[r,p] = corr(SSD2(:),dat3(:),'rows','complete','Type','Spearman');
title(['r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))])
xlabel('Change in FR to Replays of Different Arms')
ylabel('FR to Non-Local HP Representation while running')
set(gca,'FontSize',18)
%local high veloccity
subplot(1,2,2); hold on

dat2 = dat(:,:,2,1,ilab,1);
dat2(isinf(dat2)) = NaN;
% dat3 = (dat2-median(dat2,2));
% SSD2 = (SSD-median(SSD,2));
dat3 = nanzscore(dat2,[],2); % zscored
SSD2 = zscore(SSD,[],2);
plot(SSD2(:),dat3(:),'ok','MarkerSize',10,'LineWidth',3)
xl = get(gca,'xlim');
xx = SSD2(:); yy = dat3(:);
lm = polyfit(xx(~isnan(yy) & ~isnan(xx)),yy(~isnan(yy) & ~isnan(xx)),1);
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',2)
[r,p] = corr(SSD2(:),dat3(:),'rows','complete','Type','Spearman');
xlabel('Change in FR to Replays of Different Arms')
ylabel('FR to Local HP Representation while running')
title(['r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))])
set(gca,'FontSize',18)
suptitle('Across PFC cells, zscored for each cell')
set(gcf,'Position',[2126         -40        1351         725])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\AcrossDays_NonLocalVLocal_zscored_' label{ilab}])
%
%non-local high veloccity, partialling out local
figure; hold on
subplot(1,2,1); hold on
dat1 = dat(:,:,2,2,2,1);
dat2 = dat(:,:,2,1,2,1);
dat2(isinf(dat2)) = NaN;
dat4 = (dat2-median(dat2,2)); %zscored local
dat3 = (dat1-median(dat1,2)); %zscored nonlocal
SSD2 = (SSD-median(SSD,2));
exclude = isnan(dat3(:,1));
SSD2 = SSD2(~exclude,:); dat3 = dat3(~exclude,:); dat4 = dat4(~exclude,:);
% dat3 = zscore(dat1,[],2); % zscored non-local
% dat4 = zscore(dat2,[],2); % zscored local
% SSD2 = zscore(SSD,[],2);

SSD3 = SSD2(:); dat5 = dat3(:); 
dat6 = dat4(:); %partial this out

z = dat6;
z1 = [ones(size(z,1),1) z];
x = [SSD3 dat5];
resid = x - z1*(z1 \ x);
% [r,p] = corr(resid(:,1),resid(:,2));
plot(resid(:,1),resid(:,2),'ok','MarkerSize',10,'LineWidth',3)
xl = get(gca,'xlim');
xx = resid(:,1); yy = resid(:,2); xx = xx(:); yy = yy(:);
lm = polyfit(xx(~isnan(yy) & ~isnan(xx)),yy(~isnan(yy) & ~isnan(xx)),1);
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',2)
[r,p] = partialcorr(SSD3,dat5,dat6);
xlabel('Residual(FR to Replays)')
ylabel('Residual(FR to Non-Local HP while running)')
title(['Partialling out Local, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))])
set(gca,'FontSize',18)

subplot(1,2,2); hold on

% lm2 = fitlm(dat6,dat5,'linear');
% resdat6_2 = lm2.Residuals.Raw;
% lm2 = fitlm(SSD3,dat5,'linear');
% resSSD3_2 = lm2.Residuals.Raw;
% [r,p] = corr(resSSD3_2,resdat6_2);
% plot(resSSD3_2,resdat6_2,'ok','MarkerSize',10,'LineWidth',3)

z = dat5;
z1 = [ones(size(z,1),1) z];
x = [SSD3 dat6];
resid = x - z1*(z1 \ x);
% [r,p] = corr(resid(:,1),resid(:,2))
plot(resid(:,1),resid(:,2),'ok','MarkerSize',10,'LineWidth',3)
xl = get(gca,'xlim');
xx = resid(:,1); yy = resid(:,2); xx = xx(:); yy = yy(:);
lm = polyfit(xx(~isnan(yy) & ~isnan(xx)),yy(~isnan(yy) & ~isnan(xx)),1);
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',2)
[r,p] = partialcorr(SSD3,dat6,dat5);
xlabel('Residual(FR to Replays)')
ylabel('Residual(FR to Local HP while running)')
title(['Partialling out Non-Local, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))])
set(gca,'FontSize',18)
set(gcf,'Position',[2126         -40        1351         725])
suptitle('Across PFC cells, Partial Correlations (on Median subtracted data)')
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\AcrossDays_NonLocalVLocal_MedianSubtracted_Partial_' label{ilab}])

%


%non-local high veloccity, partialling out local
figure; hold on
subplot(1,2,1); hold on
dat1 = dat(:,:,2,2,2,1);
dat2 = dat(:,:,2,1,2,1);
dat2(isinf(dat2)) = NaN;
dat3 = zscore(dat1,[],2); % zscored non-local
dat4 = zscore(dat2,[],2); % zscored local
% dat4 = (dat2-median(dat2,2)); %zscored local
% dat3 = (dat1-median(dat1,2)); %zscored nonlocal
SSD2 = (SSD-median(SSD,2));
exclude = isnan(dat3(:,1));
SSD2 = SSD2(~exclude,:); dat3 = dat3(~exclude,:); dat4 = dat4(~exclude,:);

% SSD2 = zscore(SSD,[],2);

SSD3 = SSD2(:); dat5 = dat3(:); 
dat6 = dat4(:); %partial this out

z = dat6;
z1 = [ones(size(z,1),1) z];
x = [SSD3 dat5];
resid = x - z1*(z1 \ x);
% [r,p] = corr(resid(:,1),resid(:,2));
plot(resid(:,1),resid(:,2),'ok','MarkerSize',10,'LineWidth',3)
xx = resid(:,1); yy = resid(:,2); xx = xx(:); yy = yy(:);
lm = polyfit(xx(~isnan(yy) & ~isnan(xx)),yy(~isnan(yy) & ~isnan(xx)),1);
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',2)
[r,p] = partialcorr(SSD3,dat5,dat6);
xlabel('Residual(FR to Replays)')
ylabel('Residual(FR to Non-Local HP while running)')
title(['Partialling out Local, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))])
set(gca,'FontSize',18)

subplot(1,2,2); hold on
z = dat5;
z1 = [ones(size(z,1),1) z];
x = [SSD3 dat6];
resid = x - z1*(z1 \ x);
% [r,p] = corr(resid(:,1),resid(:,2))
plot(resid(:,1),resid(:,2),'ok','MarkerSize',10,'LineWidth',3)
xx = resid(:,1); yy = resid(:,2); xx = xx(:); yy = yy(:);
lm = polyfit(xx(~isnan(yy) & ~isnan(xx)),yy(~isnan(yy) & ~isnan(xx)),1);
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',2)
[r,p] = partialcorr(SSD3,dat6,dat5);
xlabel('Residual(FR to Replays)')
ylabel('Residual(FR to Local HP while running)')
title(['Partialling out Non-Local, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))])
set(gca,'FontSize',18)
set(gcf,'Position',[2126         -40        1351         725])
suptitle('Across PFC cells, Partial Correlations (on zscored data)')
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\AcrossDays_NonLocalVLocal_Zscore_Partial_' label{ilab}])
end