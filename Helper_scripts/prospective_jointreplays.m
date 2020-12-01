function prospective_jointreplays(thisdir,label,area,jointpart)
% num_shuffles,

load(thisdir,[label '_replay_stnd'],'other_cells','hp_cells','hpinterneurons','spikedata',[label '_replay_singlebothjoint'],...
    [label  '_replay_jointreplayarm'],[label  '_replay_jointjumptime'],[label  '_replay_post'],[label  '_replay_postseqindex'],...
    [label  '_replay_shuffle_p'])

eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
% eval(['Eventp = ' label '_replay_shuffle_p;'])
if jointpart==1
    eval(['Event = ' label '_replay_stnd(:,1);'])
elseif jointpart == 2
    eval(['Event = ' label '_replay_stnd(:,2);'])
end
eval(['JointArm = ' label  '_replay_jointreplayarm;'])
eval(['JointTime = ' label  '_replay_jointjumptime;'])
% eval(['EventEnd = ' label '_replay_stnd(:,2);'])
clear([label '_replay_stnd'])
if area == 1
    p = other_cells; 
elseif area == 2
    p = setdiff(hp_cells,hpinterneurons);
end
clear other_cells hp_cells hpinterneurons

% Event(singlebothjoint~=3 | Eventp>=.05) = [];
% JointArm(singlebothjoint~=3 | Eventp>=.05,:) = [];
% JointTime(singlebothjoint~=3 | Eventp>=.05,:) = [];

Event(singlebothjoint~=3) = [];
JointArm(singlebothjoint~=3,:) = [];
JointTime(singlebothjoint~=3,:) = [];


spk = spikedata(ismember(spikedata(:,2),p),1:2); 


JFR = NaN(size(Event,1),length(p));


if jointpart==1
    for iE = 1:size(Event,1)
        for icell = 1:length(p)        
%             JFR(iE,icell) = sum(spk(:,2)==p(icell) & spk(:,1)>=Event(iE) & spk(:,1)<EventEnd(iE))./(EventEnd(iE)-Event(iE)); %changed 10/22 to change in FR 
            JFR(iE,icell) = sum(spk(:,2)==p(icell) & spk(:,1)>=Event(iE) & spk(:,1)<JointTime(iE))./(JointTime(iE)-Event(iE)); %changed 10/22 to change in FR 
%         JFR(iE,icell) = (sum(spk(:,2)==p(icell) & spk(:,1)>=Event(iE) & spk(:,1)<JointTime(iE))./(JointTime(iE)-Event(iE)))- ...
%             (sum(spk(:,2)==p(icell) & spk(:,1)<Event(iE) & spk(:,1)>=(Event(iE)-.2))./.2);     
        end
    end
elseif jointpart==2
    for iE = 1:size(Event,1)
        for icell = 1:length(p)        
            JFR(iE,icell) = sum(spk(:,2)==p(icell) & spk(:,1)<=Event(iE) & spk(:,1)>JointTime(iE))./(Event(iE)-JointTime(iE));
        end
    end
end

% new
% JFR(JFR == 0) = NaN;

p1 = NaN(3,size(JFR,2));
dfr = NaN(size(JFR,2),3);
dfrs = NaN(size(JFR,2),3,num_shuffles);
cellnums = NaN(3,1);
% propsig = cellnums;
armtowards = NaN(size(JFR,2),3,2);
for iarm = 1:3
    
    if jointpart == 1
        touse = JointArm(:,1)==iarm;
        Labels = JointArm(touse,2);
    elseif jointpart == 2
        touse = JointArm(:,2)==iarm;
        Labels = JointArm(touse,1);
    end
    JFRarm = JFR(touse,:);
    aa = unique(Labels);
    hh = min(hist(Labels,aa));
    if length(aa)==2 && hh>2
        cdiff = abs(nanmean(JFRarm(Labels==aa(1),:),1)-nanmean(JFRarm(Labels==aa(2),:),1))';
        armtowards(:,iarm,1) = nanmean(JFRarm(Labels==aa(1),:),1);
        armtowards(:,iarm,2) = nanmean(JFRarm(Labels==aa(2),:),1);
%         sh = randi(length(Labels),length(Labels),num_shuffles);    
%         cdiffs = NaN(size(cdiff,1),num_shuffles);
%         for j = 1:num_shuffles
%             cdiffs(:,j) = abs(nanmean(JFRarm(Labels(sh(:,j))==aa(1),:),1)-nanmean(JFRarm(Labels(sh(:,j))==aa(2),:),1));           
%         end
        x = sum(cdiffs>=cdiff,2);
        p0 = (x+1)/(num_shuffles+1);
        cellnum = sum(sum(JFRarm,1)~=0);

        cdiff(sum(JFRarm,1)==0) = NaN;
        dfr(:,iarm) = cdiff;    
        dfrs(:,iarm,:) = cdiffs;
        p1(iarm,:) = p0; p1(iarm,sum(JFRarm,1)==0) = NaN;
        cellnums(iarm,1) = cellnum;
%         propsig(iarm,1) = sum(p0<.1)/cellnum;    
    end
end

toward1 = NaN(size(armtowards,1),3); toward2 = toward1; toward3 = NaN(size(armtowards,1),6);
 for icell = 1:size(armtowards,1)
        dat = squeeze(armtowards(icell,:,:,:));
        toward1(icell,1) = nanmean([dat(2,1)-dat(2,2) dat(3,1)-dat(3,2)]);
        toward1(icell,2) = nanmean([dat(1,1)-dat(1,2) dat(3,2)-dat(3,1)]);
        toward1(icell,3) = nanmean([dat(1,2)-dat(1,1) dat(2,2)-dat(2,1)]);
%         toward1(icell,sum(isnan(dat),2)>0) = NaN;
        
        toward2(icell,1) = nanmean([dat(2,1) dat(3,1)]);
        toward2(icell,2) = nanmean([dat(1,1) dat(3,2)]);
        toward2(icell,3) = nanmean([dat(1,2) dat(2,2)]);
        toward3(icell,:) = (dat(:)-nanmin(dat(:)))./(nanmax(dat(:))-nanmin(dat(:)));
end
% propsigall = nansum(nansum(p1<.05))/nansum(cellnums);
% meddiffarm = nanmedian(repmat(dfr,[1 1 num_shuffles])-dfrs,3);

% lap by lap and 1/2 and 1/2

% [dfr,dfrs,p1,toward1,toward2,toward3]
prospectiveJR_armcombs = toward3;

save(thisdir,'prospectiveJR_armcombs','-append');