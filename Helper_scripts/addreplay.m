function replayqual = addreplay(D,st,nd,posind,type,replayqual,arm,In,Out,times,jj,replaybins,numshuff)

% [replayseqindex,replayseq,replayseqtimes, replaymaxjumpdist,replaydistance,replaycorr,replaytimes,replayarm,jointreplayarm,,replayjointsingleboth,replaydirbias,replaypost,replaypostindex]

replayqual.currentnumber = replayqual.currentnumber+1;
replayqual.singlebothjoint = cat(1,replayqual.singlebothjoint,type);

replayqual.times = cat(1,replayqual.times,times(st:nd)');
replayqual.stnd = cat(1,replayqual.stnd,[times(st) times(nd)]);

%this has zeros
post = zeros(length(st:nd),size(D,2));
post(:,posind) = D(replaybins(st:nd),posind);
replayqual.post = cat(1,replayqual.post,post);

%this is tight post
post2 = D(replaybins(st:nd),posind);
replayqual.postseqindex = cat(1,replayqual.postseqindex,max([0 max(replayqual.postseqindex)])+length(st:nd));


[wc,p] = columncycleshuffle(post2,numshuff);
replayqual.wc = cat(1,replayqual.wc,wc);
replayqual.shuffle_p = cat(1,replayqual.shuffle_p,p);

Mat = [D(replaybins(st:nd),posind)'./((ones(size(D(replaybins(st:nd),posind),2),1)*[sum(D(replaybins(st:nd),posind),2)]'))];
I2=sum([1:size(Mat,1)]'*ones(1,size(Mat,2)).*Mat);
% M = Mat';
% [x,y] = find(M==repmat(max(M,[],2),[1 size(M,2)]));
% seqmax = NaN(size(I2));
% seqmax(x) = y;

seqmax = [1:size(Mat,1)]*(Mat./sum(Mat,1));


armcov = range(seqmax)/length(posind);
replayqual.armcov = cat(1,replayqual.armcov,armcov);
replayqual.distance = cat(1,replayqual.distance,range(seqmax));
replayqual.maxjump = cat(1,replayqual.maxjump,max(abs(diff(seqmax)))./length(posind));

% 1 is pure inbound, -1 is pure outbound
DirBias = (sum(sum(In(replaybins(st:nd),posind)))-sum(sum(Out(replaybins(st:nd),posind))))./...
    (sum(sum(In(replaybins(st:nd),posind)))+sum(sum(Out(replaybins(st:nd),posind))));
replayqual.dirbias = cat(1,replayqual.dirbias,DirBias);
clear In Out

jumpi = false(size(I2'));
if type<3
    replayqual.replayarm = cat(1,replayqual.replayarm,arm);
    replayqual.jointreplayarm = cat(1,replayqual.jointreplayarm,[NaN NaN]);        
    replayqual.jointjumptime = cat(1,replayqual.jointjumptime,jj(1));
elseif type==3
    replayqual.replayarm = cat(1,replayqual.replayarm,NaN);
    replayqual.jointreplayarm = cat(1,replayqual.jointreplayarm,arm);
    replayqual.jointjumptime = cat(1,replayqual.jointjumptime,jj(1));
    jumpi(jj(2)) = true;    
end
replayqual.seqjointjump = cat(1,replayqual.seqjointjump,jumpi);



replayqual.seqmax = cat(1,replayqual.seqmax,seqmax');
replayqual.seq = cat(1,replayqual.seq,I2');
replayqual.seqtimes = cat(1,replayqual.seqtimes,times(st:nd)');
replayqual.seqreplayindex = cat(1,replayqual.seqreplayindex,ones(length(st:nd),1)*replayqual.currentnumber);


