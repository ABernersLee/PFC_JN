% getting other days out for accuracy and learning curve

figlab = 'F:\XY_matdata\Figures\ForPaperReviews\AllCells_AllArmEvents\Behavior\';
%% make params
% params.daydir % this is the folder where raw data is
% params.Run_Times %get the run times not sleep times
% params.Rat_Name % letter of rat
% params.Date % actuall date of behavior

load('F:\XY_matdata\OtherDays_Behavior\dirs_otherdays.mat','dirs')

cd(dirs.homedir)
cd ..
thisdir = cd;
d2 = dir;
% id = 3:14

% params.daydir = [thisdir '\' d2(id).name '\'];
% params.Rat_Name = d2(id).name(1);
% params.Run_Times = []./1e6;

id = 3;
params.Track_Type = 2;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
params.Date = 20090603; %C1 bad position data, need to do manually 
% params.Run_Times = [34750453517  34824207835; 35068917355 36511653825; 
params.Run_Times = [38542606722 38951425661;  39332513392  39981980677]./1e6; % dont know which are real
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')


%%
id = 4;
params.Track_Type = 2; params.numarms = 3;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
params.Date = 20090605; %C3
params.Run_Times = [12149574077 15064519815; 17208846155 20212422200]./1e6;
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')

id = 5;
params.Track_Type = 2; params.numarms = 3;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
params.Date = 20090608; %C6
params.Run_Times = [264875132 1810570731; 2217200193 4695418310;]./1e6;
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')

id = 6;
params.Track_Type = 2; params.numarms = 3;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
params.Date = 20090716; %D4
params.Run_Times = [40012303441 41383010035; 41747591438 42900062514; 43370847661 46176238938]./1e6;
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')

id = 7;
params.Track_Type = 2; params.numarms = 3;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
params.Run_Times = [734870146 1303260156; 1752302183 1984511366; 2102636740 4249646418; 4462716794 5450010551; 5696391433 7378967307]./1e6;
params.Date = 20090717; %D5
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')

id = 8;
params.Track_Type = 2; params.numarms = 3;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
params.Run_Times = [257868106 4761251301]./1e6;
params.Date = 20100317; %E2
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')

id = 9;
params.Track_Type = 2; params.numarms = 3;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
params.Run_Times = [365946175 4867998383; ]./1e6;
params.Date = 20100320; %E5
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')

id = 10;
params.Track_Type = 2; params.numarms = 3;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
params.Run_Times = [29978097416 33832234929; 34276185455 35688457109]./1e6;
params.Date = 20100321; %E6
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')

id = 11;
params.Track_Type = 2; params.numarms = 3;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
params.Run_Times = [278490440 2613565920; 5129791535 12373855390]./1e6;
params.Date = 20100322; %E7
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')

id = 12;
params.Track_Type = 2; params.numarms = 3;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
params.Run_Times = [34832085260 36814297975; 37171241655 39573708684]./1e6;
params.Date = 20100323; %E8
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')

id = 13;
params.Track_Type = 2; params.numarms = 3;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
% params.Run_Times = [7957001024 8656806385; 8724757982 11940599827;
params.Run_Times = [12188172251 16943532406; 17460477857 19923832274]./1e6;
%     ; 20028210700 24614442955]./1e6;
params.Date = 20111219; %I2
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')

id = 14;
params.Track_Type = 2; params.numarms = 3;
params.daydir = [thisdir '\' d2(id).name '\'];
params.Rat_Name = d2(id).name(1);
% params.Run_Times = [15418499759 19020750654; 
params.Run_Times = [19377939715 22681787953; 23328978866 27720146293]./1e6;
%     27768016242 31520137130]./1e6;
params.Date = 20111220; %I3
save([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')



%% extrct data and get laps - id = 3 manually


load([dirs.homedir '\' params.Rat_Name '_' num2str(params.Date)],'params')
cd(params.daydir)
d2 = dir('*.mpg');
rawpos = []; runs = [];
for id = 1:size(d2,1)-1
    Video = VideoReader([d2(id).name]);
    thisframe = 1;
%     starttime = Video.CurrentTime;
    Nframes = Video.FrameRate*Video.Duration;

    p = [];
    while Video.CurrentTime<Video.Duration
        Frame = read(Video,thisframe);
        [x,y] = find(max(max(Frame))==Frame);
        p = cat(1,p,[Video.CurrentTime+params.Run_Times(id,1) mean(x) mean(y)]);        
        thisframe = thisframe+1;
        if mod(thisframe,1000)==0
            disp(['Done with ' num2str(thisframe) ' frames, ' num2str(round(100*(thisframe./Nframes))) '%'])
        end
    end
    rawpos = cat(1,rawpos,p);
    runs = cat(1,runs,id*ones(size(p,1),1));
end
rundat = false(size(rawpos,1),max(runs));
for irun= 1:max(runs)
    rundat(runs==irun,irun) = true;
end
rawpos(:,2)=rawpos(:,2)-(min(rawpos(:,2))-1);
rawpos(:,3)=rawpos(:,3)-(min(rawpos(:,3))-1);

params.ident = [params.Rat_Name '_' num2str(params.Date)];
save([dirs.spikedatadir '\' params.ident '.mat'],'rawpos','rundat','params','-append')

%% extrct data and get laps 


load('F:\XY_matdata\OtherDays_Behavior\dirs_otherdays.mat','dirs')
ExtractData(dirs)

ProcessData(dirs)

%% get accuracy for days without neural data

load('F:\XY_matdata\OtherDays_Behavior\dirs_otherdays.mat','dirs')
cd(dirs.homedir)
d2 = dir('*.mat');
otheracc = NaN(size(d2,1),10);
for id = 1:size(d2,1)
    armacc = get_behavior_accuracy(d2(id).name);  
    load(d2(id).name,'laps_coverspace','laps_singlepass','params')
    %reward rate, inward choices, alternation accuracy, min number of arm visits, number of passes

    if armacc(1,1)~=1
        armacc2 = [3 armacc(1,1) NaN NaN; armacc];
    else
        armacc2 = armacc;
    end
    
    ch = nchoosek(1:3,2);
    paths = [ch; ch(:,2:-1:1)];
    numonpaths = NaN(size(paths,1),1);
    for ipath = 1:size(paths,1)
        numonpaths(ipath) = sum(sum(armacc(:,1:2)==paths(ipath,:),2)==2);
    end
    
    %if only counting the first time they leave center do they return
%     inins = find(~isnan(armacc(:,4)))+1;
%     inins(inins>size(armacc,1)) = [];

    %if counting all chances to return to center
    inins = find(isnan(armacc(:,4)));
    
    armacc(:,5) = NaN;
    armacc(inins,5) = armacc(inins,3);
    

    
    otheracc(id,:) = [nanmean(armacc(:,3)) nanmean(armacc(:,5)) nanmean(armacc(:,4)) floor(nansum(armacc(:,4))/2) min(histc([armacc(:,1);armacc(end,2)],1:3)) min(numonpaths) max(laps_singlepass) double(params.Rat_Name) params.Date false];        
end
   
%% then also days with neural data

load('F:\XY_matdata\dirs.mat','dirs')
cd(dirs.homedir)
d2 = dir('*.mat');
theacc = NaN(size(d2,1),10);
for id = 1:size(d2,1)
    armacc = get_behavior_accuracy(d2(id).name);  
    load(d2(id).name,'laps_coverspace','laps_singlepass','params')
    %reward rate, inward choices, alternation accuracy, number of full trials (C,L,C,R), min number of arm visits, number of path visits, number of passes, rat, date, dataday

    %if only counting the first time they leave center do they return
%     inins = find(~isnan(armacc(:,4)))+1;
%     inins(inins>size(armacc,1)) = [];
     if armacc(1,1)~=1
        armacc2 = [3 armacc(1,1) NaN NaN; armacc2];
    else
        armacc2 = armacc;
    end    
    ch = nchoosek(1:3,2);
    paths = [ch; ch(:,2:-1:1)];
    numonpaths = NaN(size(paths,1),1);
    for ipath = 1:size(paths,1)
        numonpaths(ipath) = sum(sum(armacc(:,1:2)==paths(ipath,:),2)==2);
    end

    %if counting all chances to return to center
    inins = find(isnan(armacc(:,4)));
    
    armacc(:,5) = NaN;
    armacc(inins,5) = armacc(inins,3);
    
   
    
    theacc(id,:) = [nanmean(armacc(:,3)) nanmean(armacc(:,5)) nanmean(armacc(:,4)) floor(nansum(armacc(:,4))/2) min(histc([armacc(:,1);armacc(end,2)],1:3)) min(numonpaths) max(laps_singlepass) double(params.Rat_Name) params.Date true];    
end
%% combine
acc1 = [otheracc;theacc];
[acc2,ord] = sort(acc1(:,end-1));
acc = acc1(ord,:);

%% plot and stats of trends


rats = unique(acc(:,8));
alldays = NaN(length(rats),8,6); %rats, day, type
c = varycolor(length(rats));
typelab = {'Overall accuracy';'Inbound accuracy';'Outbound accuracy';'Completed trials (C,R,C,L)';'Minimum times visited an arm';'Mininum times doing a path (C->R)';'Number of passes'};
for itype = 1:7
    clear f
    figure; hold on
    for irat = 1:length(rats)
        f(irat,1) = plot(acc(acc(:,8)==rats(irat),itype),'.-','Color',c(irat,:),'MarkerSize',15);
        plot(find(acc(:,8)==rats(irat) & acc(:,end)==1)-(find(acc(:,8)==rats(irat),1,'first')-1), acc(acc(:,8)==rats(irat) & acc(:,end)==1,itype),'*','Color',c(irat,:),'MarkerSize',20);
        alldays(irat,1:sum(acc(:,8)==rats(irat)),itype) = acc(acc(:,8)==rats(irat),itype);
    end
    if itype <4
       xl = get(gca,'xlim');
       plot(xl,[.5 .5],'--k')
    end
    legend(f,{'C';'D';'E';'I'},'Location','best')
    xlabel('Day')
    ylabel(typelab{itype})
%     errorbar(1:8,nanmean(alldays(:,:,itype)),nanstd(alldays(:,:,itype))./sum(~isnan(alldays(:,:,itype))),'k')
end

%%

rats = unique(acc(:,8));
alldays = NaN(length(rats),8,6); %rats, day, type
c1 = varycolor(length(rats));
cc = cat(3,c1,.5*ones(size(c1)));
c = mean(cc,3);
typelab = {'Overall accuracy';'Inbound accuracy';'Outbound accuracy';'Completed trials (C,R,C,L)';'Minimum times visited an arm';'Mininum times doing a path (C->R)';'Number of passes'};
for itype = 1:7
    clear f
    fromto = NaN(length(rats),2);
    figure; hold on
    for irat = 1:length(rats)
        f(irat,1) = plot(acc(acc(:,8)==rats(irat),itype),'.-','Color',c(irat,:),'MarkerSize',15); %'Color',[.5 .5 .5],'MarkerSize',15);
        plot(find(acc(:,8)==rats(irat) & acc(:,end)==1)-(find(acc(:,8)==rats(irat),1,'first')-1), acc(acc(:,8)==rats(irat) & acc(:,end)==1,itype),'*','Color',c(irat,:),'MarkerSize',20); %'Color',[.5 .5 .5],'MarkerSize',20);
        alldays(irat,1:sum(acc(:,8)==rats(irat)),itype) = acc(acc(:,8)==rats(irat),itype);
        dat = acc(acc(:,8)==rats(irat),itype);
        fromto(irat,:) = [dat(1) dat(end)];
    end
    xlim([.5 8.5])
    if itype <4
       xl = get(gca,'xlim');
       plot(xl,[.5 .5],'--k')
       ylim([-.1 1.1])
    end
    xlabel('Day')
    
    ylabel(typelab{itype})
    plind = find(sum(~isnan(alldays(:,:,itype)))>1,1,'last');
    errorbar([1:1:plind]+.15,nanmean(alldays(:,1:plind,itype)),nanstd(alldays(:,1:plind,itype))./sum(~isnan(alldays(:,1:plind,itype))),'k','LineWidth',2)
    
    alldat = alldays(:,:,itype); alldaydat = repmat([1:8],[4 1]);
    [r,p] = corr(alldaydat(:),alldat(:),'rows','complete','type','Kendall');
    title(['R = ' num2str(round(r,4,'significant')) ', P = ' num2str(round(p,4,'significant')) ', from ' num2str(round(mean(fromto(:,1)*100),2)) ' to ' num2str(round(mean(fromto(:,2)*100),2)) '%'])
        
    legend(f,{'C';'D';'E';'I'},'Location','best')
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([figlab '_' num2str(itype)])
end

%%
rats = unique(acc(:,8));
alldays = NaN(length(rats),8,6); %rats, day, type
c = varycolor(length(rats));
typelab = {'Overall accuracy';'Inbound accuracy';'Outbound accuracy';'Completed trials (C,R,C,L)';'Minimum times visited an arm';'Mininum times doing a path (C->R)';'Number of passes'};
for itype = 1:7
    clear f
    figure; hold on
    for irat = 1:length(rats)
%         f(irat,1) = plot(acc(acc(:,8)==rats(irat),itype),'.-','Color',[.5 .5 .5],'MarkerSize',15);
        f(irat,1) = plot(find(acc(:,8)==rats(irat) & acc(:,end)==1)-(find(acc(:,8)==rats(irat),1,'first')-1), acc(acc(:,8)==rats(irat) & acc(:,end)==1,itype),'.-','Color',[.5 .5 .5],'MarkerSize',20);
        alldays(irat,find(acc(:,8)==rats(irat) & acc(:,end)==1)-(find(acc(:,8)==rats(irat),1,'first')-1),itype) = acc(acc(:,8)==rats(irat) & acc(:,end)==1,itype);
    end
    if itype <4
       xl = get(gca,'xlim');
       plot(xl,[.5 .5],'--k')
       ylim([0 1])
    end
    
    xlabel('Day')
    ylabel(typelab{itype})
    errorbar([1:8]+.1,nanmean(alldays(:,:,itype)),nanstd(alldays(:,:,itype))./sum(~isnan(alldays(:,:,itype))),'k','LineWidth',2)
    legend(f,{'C';'D';'E';'I'},'Location','best')
end




    