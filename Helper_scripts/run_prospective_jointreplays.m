cd('E:\XY_matdata\AllDays')
d2 = dir('*.mat');
num_shuffles = 80;
num_hp_shuffles = 40;
%%
% singlecell_jointreplay_stemdiff(thisdir,label,num_shuffles,area)
% predictjointreplays.m

 num_shuffles = 1000;
    for area = 1:2
        hs = NaN(3,3,11); ps = NaN(11,3); ps2 = NaN(11,3); 
        for id = 1:11    
                [ps(id,1),ps2(id,:),ps(id,2),ps(id,3),hs(:,:,id)] = make_jointreplayeventtriggeredmat(d2(id).name,label,num_shuffles,area);      
        end
        save(['E:\XY_matdata\AllDaysMat\predictjointreplays_20181012_' label '_area' num2str(area)],'hs','ps','ps2','num_shuffles')
    end
    
    %%
ilabel = 1;
if ilabel == 1; label = 'RP'; elseif ilabel==2; label = 'SD'; end

%%
% datz = [];
armrepdat = [];
for id = 1:11
%     d = fwdvsrev_pfc(d2(id).name,label);
%     datz = [datz;d];
    a = get_armrepdat(d2(id).name,label);
    armrepdat = [armrepdat;a];
end

% p = signrank(datz(:,1),datz(:,2),'tail','right')
% p = signrank(datz(:,1),datz(:,2))
armrepdatz = nanzscore(armrepdat,[],2);
%%

%  Joint Replay Data
PFCreal = [];PFCshuff = [];HPreal = [];HPshuff = [];CoeffsPFC = [];CoeffsHP = [];

ps = NaN(11,3); % pfc v shuffle, pfc v hp, hp v shuff
tic
for id = 1:11    
    [ps(id,1),ps(id,2),ps(id,3),~,PFCreal1,PFCshuff1,HPreal1,HPshuff1,CoeffsPFC1,CoeffsHP1] = make_jointreplayeventtriggeredmat(d2(id).name,label,num_shuffles,num_hp_shuffles);
    PFCreal = [PFCreal;PFCreal1];
    PFCshuff = [PFCshuff; PFCshuff1];
    HPreal = [HPreal;HPreal1];
    HPshuff = [HPshuff; HPshuff1];
    CoeffsPFC = [CoeffsPFC;CoeffsPFC1];
    CoeffsHP = [CoeffsHP;CoeffsHP1];

    t = toc;
    disp(['Day ' num2str(id) ', ps: ' num2str(ps(id,:))  ', in ' num2str(t/60) ' Minutes'])        
    tic        
end
CoeffsPFCz = nanzscore(CoeffsPFC,[],2);
%%
% Joint Replay Prediction Significance
figure; hold on; histogram(nanmean(PFCshuff),'FaceColor','k'); yl = get(gca,'ylim');plot([nanmean(PFCreal) nanmean(PFCreal)],yl,'r','LineWidth',3)
p = (sum(nanmean(PFCshuff)>=nanmean(PFCreal))+1)/(size(PFCshuff,2)+1);
title(['Joint Replay Prediction: PFC vs Shuffle, p = ' num2str(p)])

figure; hold on; histogram(nanmean(HPreal),'FaceColor','k'); yl = get(gca,'ylim');plot([nanmean(PFCreal) nanmean(PFCreal)],yl,'r','LineWidth',3)
p = (sum(nanmean(HPreal)>=nanmean(PFCreal))+1)/(num_hp_shuffles+1);
title(['Joint Replay Prediction: PFC vs HP, p = ' num2str(p)])

figure; hold on; 
histogram(nanmean(HPshuff),'FaceColor','k'); 
yl = get(gca,'ylim');
plot([nanmean(nanmean(HPreal)) nanmean(nanmean(HPreal))],yl,'r','LineWidth',3)
p = (sum(nanmean(HPshuff)>=nanmean(nanmean(HPreal)))+1)/(size(HPshuff,2)+1);
title(['Joint Replay Prediction: HP vs Shuffle, p = ' num2str(p)])
%     save(['E:\XY_matdata\AllDaysMat\predictjointreplays_20181012_' label '_area' num2str(area)],'hs','ps','ps2','num_shuffles')
%% armsig replay data
armsig2 = [];
moduarmall = [];
for id = 1:11
    load(d2(id).name, [label '_moduarm'],[label '_pSSDarm'],[label '_SSDarm'],[label '_Cand_sig_modu_include'])
    eval(['armsig = [' label '_pSSDarm ' label '_SSDarm ' label '_Cand_sig_modu_include];'])
%     eval(['moduarm = [' label '_moduarm2];'])
    eval(['moduarm = [' label '_moduarm];'])
    [~,b] = max(moduarm);
    [~,ab] = max(abs(moduarm));
    armsig(:,6:7) = [b' ab'];
    combs = nchoosek(1:3,2);
    SSDnew = zeros(size(armsig,1),1);
    for ic = 1:size(combs,1)
        SSDnew = SSDnew+((moduarm(combs(ic,1),:)-moduarm(combs(ic,2),:))'.^2);    
    end
    armsig(:,8) = SSDnew;
    armsig2 = cat(1,armsig2,armsig);
    moduarmall = [moduarmall;moduarm'];
end
armsig = armsig2;
moduarmallz = nanzscore(moduarmall,[],2);
%% Prospective Coding Data & Fields & Theta cell plots
num_shuffles = 80;
num_hp_shuffles = 40;
PFCrealA = [];PFCshuffA = [];HPrealA = [];HPshuffA = [];CoeffsPFCA = [];CoeffsHPA = [];wbeta = []; Armind2= [];
ps2 = NaN(11,3); % pfc v shuffle, pfc v hp, hp v shuff
tic
toplot = false;
for id = 1:11    
    [ps2(id,1),ps2(id,2),ps2(id,3),Armind2,PFCreal2,PFCshuff2,HPreal2,HPshuff2,CoeffsPFC2,CoeffsHP2]...
        = prospective_coding_test(d2(id).name,false,label,num_shuffles,num_hp_shuffles);
    PFCrealA = [PFCrealA;PFCreal2];
    PFCshuffA = [PFCshuffA; PFCshuff2];
    HPrealA = [HPrealA;HPreal2];
    HPshuffA = [HPshuffA; HPshuff2];
    CoeffsPFCA = [CoeffsPFCA;CoeffsPFC2];
    CoeffsHPA = [CoeffsHPA;CoeffsHP2];
    
%     wbeta = cat(1,wbeta,nanmean(CoeffsPFC2,2));
%         
%         if toplot
%             [diffperc,~,~,~,pl] = prospective_coding_fields(d2(id).name);
% 
%             [pfc31,~,pfc3b1,~,~,~] = make_Theta_Seq_forXWdat(d2(id).name);
% 
%             load(d2(id).name,'dirname','other_cells')    
%             plot_prospective_fd(pl,diffperc',CoeffsPFC2,dirname,other_cells,pfc31,pfc3b1)
%         end
    t = toc;
    disp(['Day ' num2str(id) ', ps: ' num2str(ps2(id,:))  ', in ' num2str(t/60) ' Minutes'])        
    tic 
end
CoeffsPFCAz = nanzscore(CoeffsPFCA,[],2);

% Prospective Coding Significance
figure; hold on; histogram(nanmean(PFCshuffA),'FaceColor','k'); yl = get(gca,'ylim');plot([nanmean(PFCrealA) nanmean(PFCrealA)],yl,'r','LineWidth',3)
p = (sum(nanmean(PFCshuffA)>=nanmean(PFCrealA))+1)/(size(PFCshuffA,2)+1);
title(['Prospective Coding Prediction: PFC vs Shuffle, p = ' num2str(p)])

figure; hold on; histogram(nanmean(HPrealA),num_hp_shuffles,'FaceColor','k'); yl = get(gca,'ylim');plot([nanmean(PFCrealA) nanmean(PFCrealA)],yl,'r','LineWidth',3)
p = (sum(nanmean(HPrealA)>=nanmean(PFCrealA))+1)/(num_hp_shuffles+1);
title(['Prospective Coding Prediction: PFC vs HP, p = ' num2str(p)])

figure; hold on; 
histogram(nanmean(HPshuffA),'FaceColor','k'); 
yl = get(gca,'ylim');
plot([nanmean(nanmean(HPrealA)) nanmean(nanmean(HPrealA))],yl,'r','LineWidth',3)
p = (sum(nanmean(HPshuffA)>=nanmean(nanmean(HPrealA)))+1)/(size(HPshuffA,2)+1);
title(['Prospective Coding Prediction: HP vs Shuffle, p = ' num2str(p)])
%% joint replay plots
smoothsize = 30; withmodpatch = 1; 
tic
for ilabel = 1:2
    if ilabel == 1; label = 'RP'; elseif ilabel==2; label = 'SD'; end
    for id = 1:11    
         [~,~,~,~,~,~,~,~,CoeffsPFC1,~] = make_jointreplayeventtriggeredmat(d2(id).name,label,num_shuffles,num_hp_shuffles);
         [dfr1,~,~,~,t2,~] = singlecell_jointreplay_stemdiff(d2(id).name,label,num_shuffles,1,1);
        for timetrig = 1:2
            plot_PFC_JointReplay_triggered(d2(id).name,label,withmodpatch,smoothsize,timetrig,CoeffsPFC1,dfr1,t2)
        end
        t = toc;
        disp(['JR, ' label ' Day ' num2str(id) ', in ' num2str(t/60) ' Minutes'])        
        tic 
    end
end

%% Beta of prospective coding related to the armsig
if ilabel == 1; asig = .1; elseif ilabel == 2; asign = .01; end
wbeta3 = nanmean(abs(CoeffsPFCA),2);
figure; hold on;
set(gcf,'Position',[ 2073         146        1202         744])
subplot(1,3,1), hold on
h(1) = histogram(wbeta3(armsig(:,1)>=asig),[0:.1:1],'Normalization','probability');
h(2) = histogram(wbeta3(armsig(:,1)<asig),[0:.1:1],'Normalization','probability');
yl = get(gca,'ylim');
plot([nanmedian(wbeta3(armsig(:,1)>=asig)) nanmedian(wbeta3(armsig(:,1)>=asig))],yl,'b-','LineWidth',3)
plot([nanmedian(wbeta3(armsig(:,1)<asig)) nanmedian(wbeta3(armsig(:,1)<asig))],yl,'r-','LineWidth',3)
p = ranksum(wbeta3(armsig(:,1)<asig),wbeta3(armsig(:,1)>=asig),'tail','right');
legend(h,'Not Mod',['Arm Mod p < ' num2str(asig)])
tt = text(.5,.2,['p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
xlabel('Prospective Coeff')
ylabel('Probability')
set(gca,'FontSize',18)


subplot(1,3,2)
plot(log(armsig(:,2)),wbeta3,'.k','MarkerSize',15)
[r,p] = corr(armsig(:,2),wbeta3,'rows','complete');
tt = text(-10,.9,['Lin, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';    
end
tt.FontSize = 14;
[r,p] = corr(log(armsig(:,2)),wbeta3,'rows','complete');
tt = text(-10,.85,['Log, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';    
end
tt.FontSize = 14;
xlabel('log(armSSD)')
ylabel('Prospective Coeff')
set(gca,'FontSize',18)


subplot(1,3,3)
plot(armsig(:,1),wbeta3,'.k','MarkerSize',15)
[r,p] = corr(armsig(:,1),wbeta3,'rows','complete');
tt = text(.1,1.2,['Lin, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';    
end
tt.FontSize = 14;
xlabel('p-value')
ylabel('Prospective Coeff')
set(gca,'FontSize',18)

% helper_savefig(['D:\XY_matdata\Figures\ForPaper\ProspectiveCoding\AllClassifyTest_ArmSig_wbeta'])

%% FR differences in prospective coding (splitter cells)
diffperc = []; halfperc = []; whicharm = []; 
minsig = []; meansig = []; toplotfd = []; alltoward2 = []; alltoward1 = []; alltoward3 = []; alltoward4 = [];
for id = 1:11
    [percdiff,perchalf,pSSD,wa,pl] = prospective_coding_fields(d2(id).name);
    diffperc = cat(1,diffperc,percdiff');
    halfperc = cat(1,halfperc,(perchalf(:,:,2)-perchalf(:,:,1))');
    meansig = cat(1,meansig,mean(pSSD,2));
    minsig = cat(1,minsig,min(pSSD,[],2));    
    whicharm = cat(1,whicharm,wa);
    toplotfd = cat(1,toplotfd,pl);  
    toward1 = NaN(size(pl,1),3); toward2 = toward1; toward3 = NaN(size(pl,1),6); toward4 = toward1;
    for icell = 1:size(pl,1)
        dat = squeeze(max(pl(icell,:,:,:),[],2));
        toward1(icell,1) = nanmean([dat(2,1)-dat(2,2) dat(3,1)-dat(3,2)]);
        toward1(icell,2) = nanmean([dat(1,1)-dat(1,2) dat(3,2)-dat(3,1)]);
        toward1(icell,3) = nanmean([dat(1,2)-dat(1,1) dat(2,2)-dat(2,1)]);
%         toward1(icell,sum(isnan(dat),2)>0) = NaN;
        
        toward2(icell,1) = nanmean([dat(2,1) dat(3,1)]);
        toward2(icell,2) = nanmean([dat(1,1) dat(3,2)]);
        toward2(icell,3) = nanmean([dat(1,2) dat(2,2)]);
        
        toward3(icell,:) = (dat(:)-nanmin(dat(:)))./(nanmax(dat(:))-nanmin(dat(:)));
        toward4(icell,:) = nanmean(dat,2);
    end
    alltoward1 = [alltoward1;toward1];
    alltoward2 = [alltoward2;toward2];
    alltoward3 = [alltoward3;toward3];
    alltoward4 = [alltoward4;toward4];
    disp(['Done Day ' num2str(id)])
end
alltoward1z = nanzscore(alltoward1,[],2);
alltoward2z = nanzscore(alltoward2,[],2);
alltoward3z = nanzscore(alltoward3,[],2);
alltoward4z = nanzscore(alltoward4,[],2);
%%
diffpercz = nanzscore(diffperc,[],2);
%% FR differences in joint replays
area= 1; jointpart = 1;
dfr = []; dfrs = []; p1 = []; reptoward2 = []; reptoward1 = []; reptoward3 = [];
for id = 1:11
    [dfr1,dfrs1,p11,t1,t2,t3] = singlecell_jointreplay_stemdiff(d2(id).name,label,num_shuffles,area,jointpart);
    dfr = [dfr;dfr1]; % need to know which one the cell preditcs, same deal as splitter toward1 and 2
    dfrs = [dfrs;dfrs1];
    p1 = [p1;p11'];
    reptoward1 = [reptoward1;t1];
    reptoward2 = [reptoward2;t2];
    reptoward3 = [reptoward3;t3];
    disp(['Done Day ' num2str(id)])
end
p = (sum(nanmean(nanmax(dfrs,[],2))>=nanmean(nanmax(dfr,[],2)))+1)/(size(dfrs,3)+1);
p = (sum(nanmean(nanmean(dfrs,2))>=nanmean(nanmean(dfr,2)))+1)/(size(dfrs,3)+1);
dfrz = nanzscore(dfr,[],2);
reptoward2z = nanzscore(reptoward2,[],2);
reptoward2zz = nanzscore(abs(reptoward2),[],2);
reptoward1z = nanzscore(reptoward1,[],2);
reptoward3z = nanzscore(reptoward3,[],2);

%% Theta difference
pfc3 = []; pfc6 = []; pfc3b = []; Opfc3 = []; Opfc6 = []; Opfc3b = [];
for id = 1:11
    [pfc31,pfc61,pfc3b1,Opfc31,Opfc61,Opfc3b1] = make_Theta_Seq_forXWdat(d2(id).name);
    pfc3 = [pfc3;pfc31];
    pfc6 = [pfc6;pfc61];
    pfc3b = [pfc3b;pfc3b1];
    Opfc3 = [Opfc3;Opfc31];
    Opfc6 = [Opfc6;Opfc61];
    Opfc3b = [Opfc3b;Opfc3b1];
    disp(['Done Day ' num2str(id)])
end

pfc3z = nanzscore(pfc3,[],2);
pfc3bz = nanzscore(pfc3b,[],2);
pfc6z = nanzscore(pfc6,[],2);
Opfc3z = nanzscore(Opfc3,[],2);
Opfc3bz = nanzscore(Opfc3b,[],2);
Opfc6z = nanzscore(Opfc6,[],2);

%% Theta locking
mus = [];
for id = 1:11
    mu = pfc_theta_lock_byarm(d2(id).name);
    mus = [mus;mu];
end

%% Sig Things for delta FR during joint replay (and mean FR in prospective coding)

%check that the coefficients reflect the FR difference in JR (was sig,p<10^-14, now not)
[r,p] = corr(CoeffsPFCz(:),dfrz(:),'rows','complete','type','kendall')  % the coefficients of which arm has best JR prospective on are correlated with the difference in FR in that too
[r,p] = corr(CoeffsPFC(:),dfr(:),'rows','complete')
[r,p] = corr(abs(reptoward2(:)),alltoward2(:),'rows','complete')

%the prospective coding coefficients dont make sense
[r,p] = corr(CoeffsPFCA(:),diffperc(:),'rows','complete','type','kendall')
[r,p] = corr(CoeffsPFCAz(:),diffpercz(:),'rows','complete','type','kendall')
[r,p] = corr(CoeffsPFCA(:),alltoward2(:),'rows','complete','type','kendall')

[r,p] = corr(CoeffsPFCz(:),CoeffsPFCAz(:),'rows','complete','type','kendall')

%These are all sig and point to prospective coding being a theta cycle by theta cycle phenomenon (sig, p<10^-3 - p<10^-14)
[r,p] = corr(alltoward4z(:),pfc3bz(:),'rows','complete','type','kendall') % the arm the PFC cell fires most (stem), and the arm on which it has the best theta on coding for future (stem)
[r,p] = corr(alltoward2z(:),pfc3z(:),'rows','complete','type','kendall') % The arm with the most prospective coding has the most theta prospective coding. (arm and arm) SD sig! points to theta being involved in prospective coding (r=.1354, 1.4 xp<10^-5)
[r,p] = corr(alltoward4z(:),pfc3z(:),'rows','complete','type','kendall') % This is negative sig. The arm with the highest FR (stem)is not the one with the most prospective coding on it (arm)
[r,p] = corr(alltoward3z(:),pfc6z(:),'rows','complete','type','kendall') % this is positive sig too, all combinations of FR and theta prospectiveness

[r,p] = corr(alltoward3z(:),armrepdatz(:),'rows','complete','type','kendall')
[r,p] = corr(armrepdatz(:),pfc6z(:),'rows','complete','type','kendall')

% This is sig (in change FR) and shows that the arm with highest FR is the arm that gets predicted best during JR (p ~.02 with baseline of 500ms, gone with smaller baselines)
[r,p] = corr(alltoward4z(:),reptoward1z(:),'rows','complete','type','kendall') % the arm the PFC cell fires most, and the arm it predicts the best on replay
[r,p] = corr(alltoward4z(:),reptoward2z(:),'rows','complete','type','kendall') % the arm the PFC cell fires most, and the arm it predicts the best on replay
[r,p] = corr(alltoward4z(:),reptoward2zz(:),'rows','complete','type','kendall')
% ind = armsig(:,4)>.10;
ind = armsig(:,5)==1; % & armsig(:,4)>0; %not sig
a = reptoward2z(ind,:); b = alltoward4z(ind,:); %sig sometimes
a = reptoward3z(ind,:); b = alltoward3z(ind,:); %sig neg sometimes
a = reptoward3z(ind,:); b = pfc6z(ind,:); %sig neg sometimes
a = alltoward3z(ind,:); b = pfc6z(ind,:); %sig all the time
a = alltoward3z(ind,:); b = armrepdatz(ind,:); 
a = pfc6z(ind,:); b = armrepdatz(ind,:); 
a = CoeffsPFCz(ind,:); b = CoeffsPFCAz(ind,:);
[r,p] = corr(a(:),b(:),'rows','complete','type','kendall')

[r,p] = corr(alltoward3z(:),reptoward3z(:),'rows','complete','type','kendall') %neg between combinations of all prospective FR and all JR prospective
[r,p] = corr(reptoward3z(:),pfc6z(:),'rows','complete','type','kendall') %neg trending between combinations of all prospective theta and all JR prospective

[r,p] = corr(alltoward4z(:),moduarmallz(:),'rows','complete','type','kendall')
[r,p] = corr(moduarmallz(:),reptoward2zz(:),'rows','complete','type','kendall')
[r,p] = corr(moduarmallz(:),reptoward2z(:),'rows','complete','type','kendall')
[r,p] = corr(moduarmallz(:),alltoward2z(:),'rows','complete','type','kendall')
[r,p] = corr(moduarmallz(:),alltoward1z(:),'rows','complete','type','kendall')
[r,p] = corr(moduarmallz(:),diffpercz(:),'rows','complete','type','kendall')


CoeffsPFCAz = nanzscore(CoeffsPFCA,[],2); 
[r,p] = corr(moduarmallz(:),CoeffsPFCAz(:),'rows','complete') % arm that you have the best prediction on (stem) in FR model, is not the one you fire the most to in single replays


[r,p] = corr(moduarmallz(:),pfc3z(:),'rows','complete','type','kendall')
[r,p] = corr(moduarmallz(:),alltoward2z(:),'rows','complete','type','kendall')

[r,p] = corr(pfc3bz(:),reptoward2zz(:),'rows','complete','type','kendall')

%neg
[r,p] = corr(alltoward2z(:),reptoward2zz(:),'rows','complete','type','kendall') % neg corr between abs of predicted arms (arm arm)
[r,p] = corr(reptoward3z(:),alltoward3z(:),'rows','complete','type','kendall') % the jr you best predict and the arm you best predict, all combinations
[r,p] = corr(reptoward3z(:),pfc6z(:),'rows','complete','type','kendall')
%% none sig so far
musz1 = nanzscore(mus(:,:,1),[],2);
musz2 = nanzscore(mus(:,:,2),[],2);
musz3 = nanzscore(mus(:,:,3),[],2);
[r,p] = corr(musz1(:),reptoward2z(:),'rows','complete','type','kendall')
[r,p] = corr(musz2(:),reptoward2z(:),'rows','complete','type','kendall')
[r,p] = corr(musz3(:),reptoward2z(:),'rows','complete','type','kendall')
[r,p] = corr(musz1(:),alltoward2z(:),'rows','complete','type','kendall')
[r,p] = corr(musz2(:),alltoward2z(:),'rows','complete','type','kendall')
[r,p] = corr(musz3(:),alltoward2z(:),'rows','complete','type','kendall')
[r,p] = corr(musz1(:),alltoward4z(:),'rows','complete','type','kendall')
[r,p] = corr(musz2(:),alltoward4z(:),'rows','complete','type','kendall')
[r,p] = corr(musz3(:),alltoward4z(:),'rows','complete','type','kendall')
[r,p] = corr(musz1(:),dfrz(:),'rows','complete','type','kendall')
[r,p] = corr(musz2(:),dfrz(:),'rows','complete','type','kendall')
[r,p] = corr(musz3(:),dfrz(:),'rows','complete','type','kendall')

[r,p] = corr(musz1(:),moduarmallz(:),'rows','complete','type','kendall') % trending p.06, theta mrv on way in
[r,p] = corr(musz2(:),moduarmallz(:),'rows','complete','type','kendall')
[r,p] = corr(musz3(:),moduarmallz(:),'rows','complete','type','kendall') % trending p.18, theta mrv overall

%%
[r,p] = corr(diffpercz(:),pfc3z(:),'rows','complete','type','kendall') 
[r,p] = corr(reptoward2z(:),pfc3z(:),'rows','complete','type','kendall') 
[r,p] = corr(reptoward2z(:),Opfc3z(:),'rows','complete','type','kendall') 
[r,p] = corr(alltoward2z(:),Opfc3z(:),'rows','complete','type','kendall') 
[r,p] = corr(dfrz(:),Opfc3bz(:),'rows','complete','type','kendall') 
[r,p] = corr(dfrz(:),Opfc3z(:),'rows','complete','type','kendall') 
[r,p] = corr(pfc3bz(:),Opfc3z(:),'rows','complete','type','kendall') 
[r,p] = corr(pfc3z(:),Opfc3bz(:),'rows','complete','type','kendall') 
%%
rall = NaN(size(reptoward3,1),1);
for icell = 1:size(reptoward3,1)
   rall(icell,1) = corr(reptoward3(icell,:)',alltoward3(icell,:)','rows','complete','type','kendall'); 
end
%%
figure; hold on; plot(reptoward2(:),alltoward2(:),'ko')
[r,p] = corr(reptoward2(:),alltoward2(:),'rows','complete')
[r,p] = corr(reptoward2z(:),alltoward2z(:),'rows','complete')
[r,p] = corr(nanmean(reptoward2,2),nanmean(alltoward2,2),'rows','complete')


figure; hold on; plot(reptoward3(:),alltoward3(:),'ko')
figure; hold on; plot(log(reptoward3(:)),log(alltoward3(:)),'ko'); 
plot(log(nanmean(reptoward3,2)),log(nanmean(alltoward3,2)),'r+')


[r,p] = corr(reptoward3z(:),alltoward3z(:),'rows','complete','type','kendall')
[r,p] = corr(nanmean(reptoward2,2),nanmean(alltoward2,2),'rows','complete')
% figure; hold on; plot(reptoward2z(:),alltoward2z(:),'ko')
[r,p] = corr(reptoward2z(:),alltoward2z(:),'rows','complete','type','kendall');
%%


[r,p] = corr(moduarmallz(:),pfc3z(:),'rows','complete','type','kendall') 
%%


%%
CoeffsPFCz = nanzscore(abs(CoeffsPFC),[],2);
[r,p] = corr(CoeffsPFCz(:),pfc3bz(:),'rows','complete','type','kendall') 
[r,p] = corr(dfrz(:),pfc3bz(:),'rows','complete','type','kendall') 

[r,p] = corr(CoeffsPFCz(:),alltoward4z(:),'rows','complete','type','kendall') 
%%

%%
[r,p] = corr(reptoward1z(:),pfc3z(:),'rows','complete','type','kendall') %SD & RP sig neg
[r,p] = corr(reptoward2z(:),pfc3z(:),'rows','complete','type','kendall') % SD & RP sig neg
%%
[r,p] = corr(dfrz(:),pfc3bz(:),'rows','complete','type','kendall')
[r,p] = corr(reptoward3z(:),pfc6z(:),'rows','complete','type','kendall') % SD & RP sig neg
%%
[r,p] = corr(alltoward1z(:),pfc3z(:),'rows','complete','type','kendall')

[r,p] = corr(dfrz(:),alltoward4z(:),'rows','complete','type','kendall') 

[r,p] = corr(alltoward2z(:),retoward2z(:),'rows','complete','type','kendall')

[r,p] = corr(reptoward2z(:),alltoward2z(:),'rows','complete','type','kendall')

%%
% SD none significant
[r,p] = corr(reptoward1z(:),diffpercz(:),'rows','complete','type','kendall') % the joint replay you best predict and the arm on which you have the most prospective coding
[r,p] = corr(reptoward2z(:),diffpercz(:),'rows','complete','type','kendall') % the joint replay you best predict and the arm on which you have the most prospective coding
[r,p] = corr(dfrz(:),diffpercz(:),'rows','complete','type','kendall') % the joint replay on which you have the most prospective coding and the arm on which you have the most prospective coding

[r,p] = corr(reptoward1z(:),alltoward1z(:),'rows','complete','type','kendall') % the jr you best predict and the arm you best predict
[r,p] = corr(reptoward2z(:),alltoward2z(:),'rows','complete','type','kendall') % the jr you best predict and the arm you best predict


[r,p] = corr(alltoward4z(:),reptoward1z(:),'rows','complete','type','kendall') % the arm the PFC cell fires most, and the arm it predicts the best on replay
[r,p] = corr(alltoward4z(:),reptoward2z(:),'rows','complete','type','kendall') % the arm the PFC cell fires most, and the arm it predicts the best on replay
[r,p] = corr(alltoward4z(:),dfrz(:),'rows','complete','type','kendall') % the arm the PFC cell fires most, and the arm it has the best prospective coding on in joint replays
[r,p] = corr(alltoward4z(:),moduarmallz(:),'rows','complete','type','kendall')

% the # 2 is better neg correlation
[r,p] = corr(dfrz(:),reptoward1z(:),'rows','complete') % the joint replays you have the most prospective coding, and the joint replays of hwich you predict the best (should be neg or none)
[r,p] = corr(dfrz(:),reptoward2z(:),'rows','complete') % the joint replays you have the most prospective coding, and the joint replays of hwich you predict the best (should be neg or none)
[r,p] = corr(alltoward1z(:),diffpercz(:),'rows','complete') % the arms on which you have the best prospective coding, and the arms on which you have the best prediction of (should be neg or nothing)
[r,p] = corr(alltoward2z(:),diffpercz(:),'rows','complete') % the arms on which you have the best prospective coding, and the arms on which you have the best prediction of (should be neg or nothing)
%%

[r,p] = corr(dfrz(:),moduarmallz(:),'rows','complete') % SD sig at some point, arm you are most modulated by in replay, you are also most prospective on in replay

%%
[r,p] = corr(dfrz(:),moduarmallz(:),'rows','complete')
%%
[r,p] = corr(dfrz(:),diffpercz(:),'rows','complete','type','kendall')
[r,p] = corr(dfrz(:),alltoward1z(:),'rows','complete','type','kendall')
[r,p] = corr(dfrz(:),alltoward2z(:),'rows','complete','type','kendall')
[r,p] = corr(dfrz(:),alltoward4z(:),'rows','complete','type','kendall')

[r,p] = corr(moduarmallz(:),diffpercz(:),'rows','complete','type','kendall')
[r,p] = corr(moduarmallz(:),alltoward1z(:),'rows','complete','type','kendall')
[r,p] = corr(moduarmallz(:),alltoward2z(:),'rows','complete','type','kendall')
[r,p] = corr(moduarmallz(:),alltoward4z(:),'rows','complete','type','kendall')


[r,p] = corr(dfrz(:),diffpercz(:),'rows','complete','type','kendall')
[r,p] = corr(dfrz(:),diffpercz(:),'rows','complete','type','kendall')
%%
[r,p] = corr(moduarmallz(:),diffpercz(:),'rows','complete','type','kendall')
[r,p] = corr(moduarmallz(:),reptoward2z(:),'rows','complete','type','kendall')
%%
[r,p] = corr(moduarmall(:),alltoward1(:),'rows','complete')
[r,p] = corr(moduarmall(:),alltoward2(:),'rows','complete')
[r,p] = corr(dfr(:),alltoward1(:),'rows','complete')
[r,p] = corr(dfrz(:),alltoward2z(:),'rows','complete')
%%
[r,p] = corr(moduarmallz(:),alltoward1z(:),'rows','complete')
[r,p] = corr(moduarmallz(:),alltoward2z(:),'rows','complete')
[r,p] = corr(moduarmallz(:),reptoward2z(:),'rows','complete')
[r,p] = corr(moduarmallz(:),dfrz(:),'rows','complete')

%%
[r,p] = corr(armsig(:,1),nanmax(dfr,[],2),'rows','complete')
[r,p] = corr(nanmax(diffperc,[],2),nanmax(dfr,[],2),'rows','complete')
[r,p] = corr(nanmean(diffperc,2),nanmean(dfr,2),'rows','complete')
[r,p] = corr(diffperc(:),dfr(:),'rows','complete')
%%
[r,p] = corr(armsig(:,1),nanmax(diffperc,[],2),'rows','complete')
[r,p] = corr(nanmax(diffperc,[],2),nanmax(dfr,[],2),'rows','complete')
%% sig things (SD & RP!)
[r,p] = corr(armsig(:,1),nanmax(dfr,[],2),'rows','complete') %*sig diff between the extent to which you differentiate in joint replays and that you have significant differential modulation
[r,p] = corr(armsig(:,1),log(nanmax(dfr,[],2)),'rows','complete')


figure; plot(log(nanmax(diffperc,[],2)),log(nanmax(dfr,[],2)),'.k')
[r,p] = corr(log(nanmax(diffperc,[],2)),log(nanmax(dfr,[],2)),'rows','complete')
%% FR difference in prospective coding related to armsig
if ilabel == 1; asig = .1; elseif ilabel == 2; asig = .01; end

diffdat = nanmax(diffperc,[],2);
% diffdat = nanmean(diffperc,2);
figure; hold on;
set(gcf,'Position',[ 2073         146        1202         744])
subplot(1,3,1), hold on
h(1) = histogram(diffdat(armsig(:,1)>=asig),20,'Normalization','probability');
h(2) = histogram(diffdat(armsig(:,1)<asig),20,'Normalization','probability');
yl = get(gca,'ylim');
plot([nanmedian(diffdat(armsig(:,1)>=asig)) nanmedian(diffdat(armsig(:,1)>=asig))],yl,'b-','LineWidth',3)
plot([nanmedian(diffdat(armsig(:,1)<asig)) nanmedian(diffdat(armsig(:,1)<asig))],yl,'r-','LineWidth',3)
p = ranksum(diffdat(armsig(:,1)<asig),diffdat(armsig(:,1)>=asig),'tail','right');
legend(h,'Not Mod',['Arm Mod p < ' num2str(asig)])
tt = text(10,.2,['p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
xlabel('Difference In FR between Heading Arm')
ylabel('Probability')
set(gca,'FontSize',18)


subplot(1,3,2)
plot(log(armsig(:,2)),diffdat,'.k','MarkerSize',15)
[r,p] = corr(armsig(:,2),diffdat,'rows','complete');
tt = text(-10,8,['Lin, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
[r,p] = corr(log(armsig(:,2)),diffdat,'rows','complete');
tt = text(-10,7,['Log, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
xlabel('log(armSSD)')
ylabel('FR diff')
set(gca,'FontSize',18)

subplot(1,3,3)
newdat = log(diffdat); newdat(newdat==-Inf) = NaN;
plot(armsig(:,1),newdat,'.k','MarkerSize',15)
[r,p] = corr(armsig(:,1),newdat,'rows','complete');
tt = text(.3,-1,['Log, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
xlabel('p-value')
ylabel('log(FR diff)')
set(gca,'FontSize',18)
%% Is the max heading arm the same as the max modulated arm, nope
% sig for arm going to (max mod arm) for SD armsig but NOT for RP (p = .8) (hmmm)

lab = {'Arm On';'Arm Going To'};
labm = {'Max Modulation Arm';'Max Abs(Modulation) Arm'};
for iarmnow = 1 %:2
    for maxmod = 7 %6:7
%         ind = armsig(:,5) == 1 & armsig(:,4)>0; %
        
%         ind = armsig(:,5)==1 & armsig(:,4)>0 & armsig(:,3)==1 & armsig(:,1)<.05;
        ind = true(size(armsig(:,1)));
%         b = whicharm(ind,iarmnow);
        % b = setdiff(repmat(1:3',[sum(ind) 1])',whicharm(ind,:)');
%         a = armsig(ind,maxmod); % 4 is max modulation, 5 is max abs(modulation)
%         [~,b] = max(diffperc(ind,:),[],2); %track
%         [~,b] = max(alltoward2(ind,:),[],2); %track
%         [~,a] = max(dfr(ind,:),[],2); %replay
%         [~,a] = max(CoeffsPFC(ind,:),[],2); %replay
        [~,a] = max(reptoward3(ind,:),[],2); %replay
%         [~,a] = max(reptoward2(ind,:),[],2); %replay
%         [~,b] = max(alltoward2(ind,:),[],2); %track       
        [~,b] = max(alltoward3(ind,:),[],2); %track   
%         [~,b] = max(pfc3(ind,:),[],2);
        
        j = confusionmat(a,b);
%         real = (sum(j(eye(3)==1)))./(sum(sum(j)));
        real = (sum(j(eye(6)==1)))./(sum(sum(j)));

        nS = 5000;
        fake = NaN(nS,1);
        for in = 1:nS 
            j = confusionmat(a(randperm(length(a))),b); 
%             j = confusionmat(randi(6,size(a,1),1),randi(6,size(a,1),1));
%             fake(in,1) = (sum(j(eye(3)==1)))./(sum(sum(j)));
            fake(in,1) = (sum(j(eye(6)==1)))./(sum(sum(j)));
        end
        p = (sum(fake>=real)+1)/(nS+1);

        figure; histogram(fake,'FaceColor','k')
        yl = get(gca,'ylim'); hold on;
        plot([real real],yl,'r','LineWidth',3)
        title([labm{maxmod-5} ' ' lab{iarmnow} ', p = ' num2str(p)])
    end
end

%% correlations

[r,p] = corr(armsig(:,2),nanmean(abs(CoeffsPFC),2),'rows','complete')
[r,p] = corr(armsig(:,2),nanmean(abs(CoeffsPFCA),2),'rows','complete')

[r,p] = corr(armsig(:,2),nanmax(abs(CoeffsPFC),[],2),'rows','complete') %*

[r,p] = corr(armsig(:,2),nanmax(abs(CoeffsPFCA),[],2),'rows','complete')

[r,p] = corr(armsig(:,1),nanmax(abs(CoeffsPFC),[],2),'rows','complete') 
[r,p] = corr(armsig(:,1),nanmax(abs(CoeffsPFCA),[],2),'rows','complete') 
[r,p] = corr(nanmax(abs(CoeffsPFC),[],2),nanmax(abs(CoeffsPFCA),[],2),'rows','complete') 
partialcorr(nanmax(abs(CoeffsPFC),[],2),nanmax(abs(CoeffsPFCA),[],2),armsig(:,2),'rows','complete')

% does the extent to which you are prospective on arms predict prospective in replays?

[r1,p1] = corr(nanmax(abs(CoeffsPFC),[],2),nanmean(abs(CoeffsPFCA),2),'rows','complete')
[r1,p1] = corr(abs(CoeffsPFC(:)),abs(CoeffsPFCA(:)),'rows','complete')
[r1,p1] = corr(abs(CoeffsHP(:)),abs(CoeffsHPA(:)),'rows','complete')
[r1,p1] = corr(nanmax(abs(CoeffsPFC),[],2),diffdat,'rows','complete')
[r1,p1] = corr(nanmean(abs(CoeffsHP),2),nanmean(abs(CoeffsHPA),2),'rows','complete')

% across days, what is the relationship between HP and PFC prospective coding?

[r3,p3] = corr(nanmean(HPreal,2),PFCreal,'rows','complete')
[r4,p4] = corr(nanmean(HPrealA,2),PFCrealA,'rows','complete')

%% JR predict behavior?
guesses = [];
getg = NaN(11,2);
for id = 1:11
    g = JR_predict_behavior(d2(id).name,label);    
    m = floor(size(g,1)/2);
%     disp([num2str(id) ' ' num2str(sum(g(:,end))./size(g,1))])
    disp([num2str(nanmean(g(:,8))) ' '  num2str(nanmean(g(:,9)))])
%     disp([num2str(nanmean(g(g(:,5)==0,10))) ' '  num2str(nanmean(g(g(:,5)==1,10)))])
    %', pt1: ' num2str(sum(g(1:m,end))./size(g(1:m,:),1))...
     %   ', pt2: ' num2str(sum(g(end-m+1:end,end))./size(g(end-m+1:end,:),1))])
     getg(id,:) = [nanmean(g(g(:,5)==0,10)) nanmean(g(g(:,5)==1,10))];
    guesses = [guesses;g];
end
disp('total:')
disp([num2str(nanmean(guesses(guesses(:,5)==0,10))) ' '  num2str(nanmean(guesses(guesses(:,5)==1,10)))])
p = signrank(getg(:,1),getg(:,2))
signrank(guesses(:,8),guesses(:,9))
nanmean(guesses(:,8)./(guesses(:,8)+guesses(:,9)))
%%

p2 = [];
for id = 1:11
    p = pfc_replay_anovan_arm_pos(d2(id).name,label);
    p2 = [p2;p];
end