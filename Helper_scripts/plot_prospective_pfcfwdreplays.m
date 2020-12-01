function plot_prospective_pfcfwdreplays(thisdir,label,smoothsize)
load(thisdir,[label '_replay_singlebothjoint'],[label '_Cand_sig_modu_include'], ...
    [label '_PFCreplayspikes_binned'],[label '_PFCreplayspikes_list'],[label  '_replay_shuffle_p']...
    ,[label  '_replay_replayarm'],[label '_replay_stnd'],'pos','armpos',...
    [label  '_replay_corr'],[label  '_replay_dirbias'],'other_cells')


eval(['PFCreplayspikes_binned = ' label '_PFCreplayspikes_binned;'])
eval(['PFCreplayspikes_list = ' label '_PFCreplayspikes_list;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['Cand_sig_modu_include = ' label '_Cand_sig_modu_include;'])
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['Event = ' label '_replay_stnd(:,1);'])
eval(['dirbias = ' label '_replay_dirbias;'])
eval(['corrs = ' label '_replay_corr;'])

eval(['shuffle_p = ' label  '_replay_shuffle_p;'])
clear([label '_PFCreplayspikes_list'],[label '_PFCreplayspikes_binned'],[label '_Cand_sig_modu_include'])

pfc = other_cells; clear other_cells


% 1 is pure inbound, -1 is pure outbound
% dirbias = dirbias(touse);

% outbound is positive, inbound is negative
% corrs = corrs(touse);

fwd = (dirbias>0 & corrs<0) | (dirbias<0 & corrs>0);
% rev = (dirbias>0 & corrs>0) | (dirbias<0 & corrs<0);
fwd = true(size(fwd));


touse = singlebothjoint~=3 & shuffle_p<.05 & fwd;

Event = Event(touse,1);
replayarm = replayarm(touse);

[~,~,i] = histcounts(Event,pos(:,1));
% Event(i==0) = [];
replayarm(i==0) = [];
armon = armpos(i(i~=0));

sigma = 1.5;
binsize = .02;
window = [-2 2];
ind = [window(1):binsize:window(2)];

modind = [ind(ind>=0 & ind<=.2)];
baseind = [ind(ind>=-.5 & ind<=-.1)];

cs = varycolor(6);
cs(cs==1)= .7;
cs = cs([1 4 3 5 2 6],:);
armlab = {'Left';'Center';'Right'};

plotwindow = [-.6 .6];        
for icell = 1:length(pfc)
    figure; hold on
%     yl2 = NaN(3,2);
    for iarm = 1:3
        clear a
        otharms = setdiff(1:3,iarm);
        subplot(1,3,iarm), hold on
        for ioth = 1:2
            repind = armon==iarm & replayarm==otharms(ioth);
            db = PFCreplayspikes_binned(icell,:,repind);
            dat1 = squeeze(mean(db(1,:,:),3))./binsize;
            sem1 = squeeze(std(db(1,:,:)./binsize,[],3))./sqrt(size(db,3));        
        
            sz = length(dat1); x = linspace(-sz / 2, sz / 2, sz);
            gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2)); gaussFilter = gaussFilter / sum (gaussFilter);
            dat = conv(dat1',gaussFilter,'same')';
            sem = conv(sem1',gaussFilter,'same')';
            rev = dat+sem;
            fwd = dat-sem;        
                    
            p = patch([ind ind(end:-1:1)],[fwd rev(end:-1:1)],'black');
            p.FaceColor=cs((iarm-1)*2+ioth,:);
            p.EdgeColor= cs((iarm-1)*2+ioth,:);
            p.FaceAlpha=.1;
            p.EdgeAlpha=0;
            a(ioth) = plot(ind,dat,'Color',cs((iarm-1)*2+ioth,:),'LineWidth',2);
                        
        end
        axis tight
        xlim(plotwindow)               
        set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
         xlabel(['Time From Start of Replay (ms)'])
        ylabel(['FR on ' armlab{iarm} ' Arm'])    
        
%         ylim(yl) 
%         legend([a(1) a(2)],{['Forward Replay of ' armlab{otharms(1)} ' Arm'],['Forward Replay of ' armlab{otharms(2)} ' Arm']},'Location','eastoutside')       
        set(gca,'FontSize',18,'FontName','Helvetica')
%         yl2(iarm,:) = get(gca,'ylim');        
    end
    
    for iarm = 1:3
       subplot(1,3,iarm)
%        set(gca,'ylim',[min(yl2(:,1)) max(yl2(:,2))])
        yl = get(gca,'ylim');
        patch([modind modind(end:-1:1)],[ones(size(modind))*yl(1) ones(size(modind))*yl(1)+((yl(2)-yl(1))*.02)],'red','FaceAlpha',.3,'EdgeAlpha',0)
        patch([baseind baseind(end:-1:1)],[ones(size(baseind))*yl(1) ones(size(baseind))*yl(1)+((yl(2)-yl(1))*.02)],'black','FaceAlpha',.3,'EdgeAlpha',0)
    end
        
    suptitle([label ' Cell ' num2str(pfc(icell)) ', Rat ' num2str(thisdir(1)) ' ' thisdir(3:end-4)])
    set(gcf,'Position',[ 3         181        2592         779])
%     helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmReplayTriggered\FwdReplayInPosition_' label '_' thisdir(1:end-4) '_smooth' num2str(smoothsize) '_Cell' num2str(pfc(icell))])

     if ~isfolder(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\Examples\' num2str(thisdir(1:end-4)) '_Cell' num2str(pfc(icell)) '\'])
        mkdir(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\Examples\' num2str(thisdir(1:end-4)) '_Cell' num2str(pfc(icell)) '\'])
    end
   helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\Examples\' num2str(thisdir(1:end-4)) '_Cell' num2str(pfc(icell)) '\ReplayInPosition_' label])
end