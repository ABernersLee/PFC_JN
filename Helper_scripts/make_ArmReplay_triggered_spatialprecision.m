function make_ArmReplay_triggered_spatialprecision(thisdir,label)
disp('Starting make_ArmReplay_triggered_spatialprecision')

load(thisdir,[label '_replay_singlebothjoint'],[label '_replay_replayarm'], ...
   'hp_cells','hpinterneurons','other_cells',[label '_replay_seqmean'],[label '_replay_seqtimes'],[label '_replay_stnd'],...
    [label '_Cand_sig_modu_include'],[label '_pSSDarm'],[label '_replay_postseqindex'],'spikedata','armposindex','InFR','OutFR')


eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['seqmean = ' label '_replay_seqmean;'])
eval(['seqtimes = ' label '_replay_seqtimes;'])
eval(['stnd = ' label '_replay_stnd;'])
eval(['postseqindex = ' label '_replay_postseqindex;'])
eval(['p_SSD = ' label '_pSSDarm;'])

clear([label '_replay_singlebothjoint'],[label '_replay_stnd'],[label '_replay_replayarm'], ...
   [label '_replay_seqmean'],[label '_replay_seqtimes'],...
    [label '_Cand_sig_modu_include'],[label '_pSSDarm'],[label '_replay_postseqindex'])


replayarm(singlebothjoint==3) = NaN;
clear singlebothjoint

% [~,seq] = histc(seqmean,1:size(armposindex,1));

% spks_hp = spikedata(ismember(spikedata(:,2),hp),:);
% spks_pfc = spikedata(ismember(spikedata(:,2),pfc),:);
% postseqindex = [postseqindex;length(seq)+1];


% hp = hp_cells(~ismember(hp_cells,hpinterneurons)); 
% toadd = [0 find(armposindex(:,1),1,'last') find(armposindex(:,2),1,'last')];

% rep = [];
% for ii = 1:length(replayarm)
%     if ~isnan(replayarm(ii))
%         repind = postseqindex(ii):(postseqindex(ii+1)-1);
%         reppos = seq(repind)+toadd(replayarm(ii));
% %         if ~isempty(rep); if sum(seqtimes(repind(1))==rep(:,1))>0; disp('problem'); end; end
%         thesetimes = seqtimes(repind);
%         rep = [rep;thesetimes reppos];
%     end    
% end
% not using rep anymore, except for indexing them, doing them the same way
% now


% %% this way decodes all hp dim at once (too much memory, doesn't even
% come out that much faster)
% cellstouse = NaN(length(hp)-1,length(hp));
% for ihp = 1:length(hp)
%     cellstouse(:,ihp) = setdiff(hp,hp(ihp));
% end
% 
% Mat = decode_anew2(cellstouse,spikedata,InFR,OutFR,stnd(~isnan(replayarm),:));
% 
% 
% newrepos = NaN(size(Mat,2),1);
% for iarm = 1:3
%     armarmind = replayarm2==iarm;
%     armind = armposindex(:,iarm);
%     Mat2 = (Mat(armind,armarmind)./(ones(size(Mat(armind,armarmind),1),1)*sum(Mat(armind,armarmind),1)));
%     I2=sum([1:size(Mat2,1)]'*ones(1,size(Mat2,2)).*Mat2);
%     I3 = I2+find(armind,1,'first')-1;
%     [~,m] = histc(I3,1:size(armind,1));
%     newrepos(armarmind) = m;
% end    
% toc
%%
hp = hp_cells(~ismember(hp_cells,hpinterneurons)); 
stnd2 = stnd(~isnan(replayarm),:);
replayarm2 = replayarm(~isnan(replayarm));
rphp1 = [];

for ihp = 1:length(hp)
    cellstouse = setdiff(hp,hp(ihp));
    
    if ihp == 1
        [Mat,Index,Cand] = decode_anew(cellstouse,spikedata,InFR,OutFR,stnd2);
        
        aa = NaN(size(Mat,2),1); times2 = aa;
        for ievent = 1:size(replayarm2)
            aa(Index(ievent)*4+1:4*Index(ievent+1)-3) = replayarm2(ievent);
            ts = Cand(ievent,1)+.01:.005:Cand(ievent,2)+.01; % take the middle of the bin
            times2(Index(ievent)*4+1:4*Index(ievent+1)-3) = ts(1:length(Index(ievent)*4+1:4*Index(ievent+1)-3));
        end     
    else
        Mat = decode_anew(cellstouse,spikedata,InFR,OutFR,stnd2);
    end
    newrepos = NaN(size(Mat,2),1);
    for iarm = 1:3        
        armarmind = aa==iarm;
        armind = armposindex(:,iarm);
        Mat2 = (Mat(armind,armarmind)./(ones(size(Mat(armind,armarmind),1),1)*sum(Mat(armind,armarmind),1)));
        I2=sum([1:size(Mat2,1)]'*ones(1,size(Mat2,2)).*Mat2);
        I3 = I2+find(armind,1,'first')-1;
        [~,m] = histc(I3,1:size(armind,1));
        newrepos(armarmind) = m;
    end        
    rphp1 = cat(2,rphp1,newrepos);
    disp(['Done hp cell ' num2str(ihp) ' of ' num2str(length(hp))])
end
%%
cellstouse = hp;
[Mat] = decode_anew(cellstouse,spikedata,InFR,OutFR,stnd2);
rep2 = NaN(size(Mat,2),1);
for iarm = 1:3        
    armarmind = aa==iarm;
    armind = armposindex(:,iarm);
    Mat2 = (Mat(armind,armarmind)./(ones(size(Mat(armind,armarmind),1),1)*sum(Mat(armind,armarmind),1)));
    I2=sum([1:size(Mat2,1)]'*ones(1,size(Mat2,2)).*Mat2);
    I3 = I2+find(armind,1,'first')-1;
    [~,m] = histc(I3,1:size(armind,1));
    rep2(armarmind) = m;
end        
rep2 = [times2 rep2];
rep2(isnan(times2),:) = [];
%%
CandCellsSpike = [Cand(:,1)+.01 Cand(:,2)-.005];
sp = zeros(size(armposindex,1),max(spikedata(:,2)));
prec1 = NaN(size(armposindex,1),max(spikedata(:,2)));
[rep4,repord] = sortrows(rep2,1);
ra = histc(rep4(:,2),1:size(armposindex,1));
rep31 = rphp1;rep31(isnan(times2),:) = []; %new, fixed
rep3 = rep31(repord,:);
spks = spikedata;
[~,I]=histc(spks(:,1),sortrows(CandCellsSpike(:)));
spks(~mod(I,2),:)=[]; % only takes spikes within the candidate events plus .005, so should be only the 5ms after the middle of the decoded 20ms bin
spks(spks(:,1)>max(times2) | spks(:,1)<min(times2),:) = [];
 for icell = 1:max(spks(:,2))
     dat = spks(spks(:,2)==icell,1);
     [~,c] = histc(dat(:,1),rep4(:,1));
     if ismember(icell,hp)         
         p = rep3(c,icell==hp); % get the decoded positions from all cells but self 
         sp(:,icell) = histc(p,1:size(armposindex,1));  
         ra2 = histc(rep3(:,icell==hp),1:size(armposindex,1));
         prec1(:,icell) = sp(:,icell)./ra2;
     else         
         p = rep4(c,2);     
         sp(:,icell) = histc(p,1:size(armposindex,1));  
         prec1(:,icell) = sp(:,icell)./ra;
     end
      
 end
%  prec1 = sp./ra;
%%

prec = NaN(size(sp,2),2);
arms = NaN(size(sp,2),2);
for icell = 1:max(spks(:,2))
    mm = NaN(3,2,2);
    for iarm = 1:3
        dat = prec1(armposindex(:,iarm),icell);
        if rem(length(dat),3)==0
            smdat = nanmean([dat(1:3:end) dat(2:3:end) dat(3:3:end)],2);
        else
            smdat = nanmean([dat(1:3:end-rem(length(dat),3)) dat(2:3:end-rem(length(dat),3)-1) dat(3:3:end)],2);
        end
        mm(iarm,1,1) = nanmean(dat);
        mm(iarm,2,1) = nanmax(dat);
        
        mm(iarm,1,2) = nanmean(smdat);
        mm(iarm,2,2) = nanmax(smdat);        
    end
    for im = 1:2
        [~,armtouse] = max(mm(:,1,im));
        prec(icell,im) = mm(armtouse,2,im)./mm(armtouse,1,im);
    end
    [~,aaa] = max(mm(:,1,:));
    arms(icell,:) = aaa;
end

%%
% toplotind = NaN(size(prec1,2),2);
% toplotind(hp,1) = 1;
% toplotind(pfc,1) = 2;
% toplotind(pfc,2) = p_SSD(other_cells_touse(:,igroup));
% toplotind(:,3:4) = arms;
spatialprecision.prec = prec;
spatialprecision.toplot = prec1;
spatialprecision.toplotind = arms;
save(thisdir,'spatialprecision','-append')
disp('Done with make_ArmReplay_triggered_spatialprecision')
