function Fig2_Examples(dirs,savefolder)

%To Do: 
% -get ripples out and plot
% -get hp place cells out and plot spikes
% - do decoding over this period
%%
cd(dirs.homedir)
d2 = dir('*.mat');
id = 10;
thisdir = d2(id).name;
label = 'RP';

load(thisdir,[label '_replay_shuffle_p'])
load(thisdir,[label '_replay_stnd'])
load(thisdir,'binpos','pos')
load(thisdir,[label '_replay_post'],[label '_replay_postseqindex'])
load(thisdir,[label '_replay_corr'])
load(thisdir,[label '_replay_singlebothjoint'],[label '_replay_replayarm'])
eval(['ts = ' label '_replay_stnd;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['postind = ' label '_replay_postseqindex;'])
eval(['post = ' label '_replay_post;'])
eval(['wc = ' label '_replay_corr;'])
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['shuffle_p = ' label  '_replay_shuffle_p;'])
lookat = find(abs(wc)>.5 & singlebothjoint~=3 & shuffle_p<.05 & ts(:,1)>1564 & ts(:,1)<1821);
% size(lookat);
% for ilook = 1:size(lookat,1); figure; hold on; imagesc(post(postind(lookat(ilook)):postind(lookat(ilook)+1)-1,:)'); axis xy; colormap hot; title(num2str(ilook)); end
%
% figure; hold on; plot(ts(abs(wc)>.7 & replayarm==1 & shuffle_p<.001),1,'*g'); plot(ts(abs(wc)>.7 & replayarm==2 & shuffle_p<.001),1,'*r'); plot(ts(abs(wc)>.7 & replayarm==3 & shuffle_p<.001),1,'*b'); plot(pos(:,1),binpos(:,1),'k.')

tosave = lookat([1 2 18]);

% change to look like XYWs
load(thisdir,'armposindex')
postsave = post;
post1 = post(:,armposindex(:,3));
post1 = post1(:,end:-1:1);
post = [post1 post(:,armposindex(:,2)) post(:,armposindex(:,1))];

for ilook = 1:size(tosave,1)
    figure; hold on; 
%     imagesc(post(postind(tosave(ilook)):postind(tosave(ilook)+1)-1,:)'); 
    imagesc(-post(postind(tosave(ilook)):postind(tosave(ilook)+1)-1,:)'); 
    axis xy; 
    colormap gray;
    axis tight; 
    set(gca,'clim',[-.1 0])
    xt = [0:10:range(postind(tosave(ilook)):postind(tosave(ilook)+1)-1)];
    set(gca,'xtick',xt,'xticklabel',xt*5)    
    xlabel('ms')
    ylabel('Position Bins')
    set(gca,'FontSize',18,'FontName','Helvetica')
    title(num2str(ilook))
    if ~isfolder([savefolder '\Figure2\'])
        mkdir([savefolder '\Figure2\'])
    end
    helper_saveandclosefig([savefolder '\Figure2\Example_id10_Fig1_ExampleReplay' num2str(ilook) '_newXYW'])
end

%%

figure; hold on
a = subplot(3,1,1);
hold on
binpossave = binpos;
binpos1 = binpos;
binpos1(ismember(binpos,find(armposindex(:,1)))) = binpos(ismember(binpos,find(armposindex(:,1))))+sum(armposindex(:,2))+sum(armposindex(:,3))+1;
binpos1(ismember(binpos,find(armposindex(:,2)))) = binpos(ismember(binpos,find(armposindex(:,2))))-sum(armposindex(:,1))+sum(armposindex(:,3));
binpos1(ismember(binpos,find(armposindex(:,3)))) = -(binpos(ismember(binpos,find(armposindex(:,3))))-(sum(armposindex(:,2))+sum(armposindex(:,3))))+(sum(armposindex(:,3))*2);

plot(pos(:,1),binpos1(:,1),'.','MarkerSize',15,'Color','k') %,[.5 .5 .5])
% xl = [1433.76737115852,1895.17809924617];
xl = [1445,1725];
[~,~,i] = histcounts(ts(tosave,1),pos(:,1));
i(2) = i(2)+200;
plot(xl,[sum(armposindex(:,3)) sum(armposindex(:,3))],'k--')
plot(xl,[sum(armposindex(:,3))+sum(armposindex(:,2)) sum(armposindex(:,3))+sum(armposindex(:,2))],'k--')
hold on;
ic = [255 130 0;34 139 34;75 0 130]/255;
for ilook = 1:size(tosave,1)    
    plot(pos(i(ilook),1),binpos1(i(ilook),1),'d','MarkerSize',25,'Color',ic(ilook,:),'MarkerFaceColor',ic(ilook,:))
end
% axis tight
set(gca,'xlim',xl)
xt = get(gca,'xtick');
set(gca,'xticklabel',round(xt/60))
set(gca,'ylim',[1 max(binpos1)])
set(gca,'ytick',max(binpos1),'yticklabel',3.23)
ylabel('Position (m)')
% xlabel('Time (Min)')
set(gca,'FontSize',18)




load(thisdir,'spikedata','hp_cells','hpinterneurons','InFR','OutFR')
touse = hp_cells(~ismember(hp_cells,hpinterneurons));
FR = InFR(touse,:)+OutFR(touse,:);
FRsave = FR;
FR1 = FR(:,armposindex(:,3));
FR1 = FR1(:,end:-1:1);
FR = [FR1 FR(:,armposindex(:,2)) FR(:,armposindex(:,1))];

[~,m] = max(FR,[],2);
[~,mm] = sort(m);
s = spikedata(ismember(spikedata(:,2),touse),:);
m2 = NaN(max(s(:,2)),1);
m2(touse) = m;
s = spikedata(ismember(spikedata(:,2),touse) & spikedata(:,1)>=xl(1) & spikedata(:,1)<=xl(2),1:2);

b = subplot(3,1,2);
plot(s(:,1),m2(s(:,2)),'.k')
set(gca,'xlim',xl)
axis tight
xt = get(gca,'xtick');
set(gca,'xticklabel',round(xt/60))
ylabel('HP Neuron')
% xlabel('Time (min)')
set(gca,'FontSize',18)


velthresh = 5;
EstBin=0.2;
Cell_Number = touse;
Start_Time=xl(1);
End_Time=xl(2);
Cand = xl; B=cumsum([0 diff(Cand')]);
Spike = [s(:,2) s(:,1)];
TimeBins=round((End_Time-Start_Time)/EstBin)*4;
% Cell_Number=size(OutFR,1);

% Bayesian Decoding - 5ms moving step, 20ms window estimate
binspike=zeros(length(Cell_Number),TimeBins);
for i=1:4
    for CellID=1:length(Cell_Number)   
        c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
        binspike(CellID,i:4:TimeBins)=c(1:TimeBins/4);
%         binspike(CellID,i:4:TimeBins)=i*10*ones(length(1:TimeBins/4),1);         %testing
    end                      
end

FR(sum(FR~=0,2)==0,:) = 1;
FR(FR==0) = .0001;


Number_Of_Bins=size(FR,2);
term2=zeros(Number_Of_Bins,TimeBins);

if Number_Of_Bins>TimeBins    
    for TBin=1:TimeBins        
        term2(:,TBin)=prod(FR.^binspike(:,TBin),1);
    end
else
    for PosBin=1:Number_Of_Bins           
        term2(PosBin,:)=prod((FR(:,PosBin)*ones(1,TimeBins)).^binspike,1);
    end
end

term3=exp(-EstBin*sum(FR,1)')*ones(1,TimeBins);

Matrix=term2.*term3;


% Index=round(B/(EstBin));

%normalize for each time bin
Matrix = Matrix./(ones(size(Matrix,1),1)*sum(Matrix,1));

%make times when there is zero spiking equal to zeros/NaNs all the way down
Matrix(:,sum(binspike>0)==0) = 0;

% Matrixsave = Matrix;
% Matrix1 = Matrix(armposindex(:,3),:);
% Matrix1 = Matrix1(end:-1:1,:);
% Matrix = [Matrix1;Matrix(armposindex(:,2),:);Matrix(armposindex(:,1),:)];

c = subplot(3,1,3); hold on
imagesc(Matrix); 
axis xy
colormap gray
cm = get(gca,'colormap');
cm2 = cm(end:-1:1,:);
set(gca,'colormap',cm2)

set(gca,'clim',[0 .1])
% linkaxes([a b c],'x')
% xt2 = get(gca,'xtick');
% set(gca,'xtick',[])
set(gca,'xtick',max(binpos),'xticklabel',3.23)
ylabel('Position (m)')
axis tight
xl2 = get(gca,'xlim');
% plot(xl2,[find(armposindex(:,1),1,'last') find(armposindex(:,1),1,'last')],'k--')
% plot(xl2,[find(armposindex(:,2),1,'last') find(armposindex(:,2),1,'last')],'k--')
plot(xl2,[sum(armposindex(:,3)) sum(armposindex(:,3))],'k--')
plot(xl2,[sum(armposindex(:,3))+sum(armposindex(:,2)) sum(armposindex(:,3))+sum(armposindex(:,2))],'k--')
xlabel('Time (Min)')
set(gca,'FontSize',18)
set(gca,'ytick',max(binpos),'yticklabel',3.23,'xtick',[])
% colorbar('Location','westoutside')

set(gcf,'Position',[ 2124          50        1082         844])
helper_saveandclosefig([savefolder '\Figure2\Example_id10_ExampleReplayPosSpikesDecoding_XYW'])


%%
figure; hold on
a = imagesc(Matrix); 
colormap gray
cm = get(gca,'colormap');
cm2 = cm(end:-1:1,:);
set(gca,'colormap',cm2)
set(gca,'clim',[0 .1])
colorbar('Location','westoutside','Ticks',[0 .1],...
         'TickLabels',{'0','>=.1'})
     colorbar('Location','southoutside','Ticks',[0 .1],...
         'TickLabels',{'0','>=.1'})
set(gca,'FontSize',18)
set(gcf,'Position',[2220         343         348         308])    
    
helper_saveandclosefig([savefolder '\Figure2\Example_id10_ExampleReplayPosSpikesDecoding_legend'])