function make_hoverjump(thisdir,label)

load(thisdir,[label  '_replay_seqmax'],[label  '_replay_seqreplayindex'],[label  '_replay_seqtimes'])

eval(['seq = ' label  '_replay_seqmax;'])
eval(['seqtimes = ' label  '_replay_seqtimes;'])
eval(['seqindex = ' label  '_replay_seqreplayindex;'])
clear([label  '_replay_seqreplayindex'])
clear([label  '_replay_seqmax'])
CandStepSize = NaN(max(seqindex),4);
hover = []; move = []; 
seqisjump_bin = []; seqisjump_btwn = []; seqstepst_btwn = []; 
seqtimes_btwn = []; seqtimes_bin = []; seqindex_bin = []; seqindex_btwn = [];
for SeqNo = 1:max(seqindex)
    Loc=seq(seqindex==SeqNo)';
    DiffLoc=abs(diff(Loc));                   
    Mark=DiffLoc>1;
    seqisjump_bin = [seqisjump_bin; 0; Mark']; % 1 is jump, 0 is hover
    seqisjump_btwn = [seqisjump_btwn;Mark'];
    
    Check=diff([0 Mark 0]);
    seqstepst_btwn = [seqstepst_btwn; Check(1:end-1)']; %1 is start to jump, -1 is start to hover    
    jnk1 = diff(seqtimes(seqindex==SeqNo))/2;
    jnk = seqtimes(seqindex==SeqNo)+[jnk1; jnk1(end)];
    seqtimes_btwn = [seqtimes_btwn; jnk(1:end-1)];
    seqtimes_bin = [seqtimes_bin; seqtimes(seqindex==SeqNo)];
    seqindex_bin = [seqindex_bin; ones(size(Loc,2),1)*SeqNo];    
    seqindex_btwn = [seqindex_btwn; ones(size(Loc,2)-1,1)*SeqNo];    
    SS=find(Check==1);
    EE=find(Check==-1);
    if isempty(SS)
        move1=[];
        hover1=[length(Mark) sum(DiffLoc) range(Loc) SeqNo];                
    else
        move1=zeros(length(SS),4); 
        hover1=zeros(length(SS),4);
        move1(:,1)=EE-SS;
        hover1(:,1)=[SS(2:end) length(Mark)+1 ]-EE;
        move1(:,4)=SeqNo;
        hover1(:,4)=SeqNo;
        for i=1:length(SS)
            move1(i,2)=sum(DiffLoc(SS(i):EE(i)-1));
            move1(i,3)=range(Loc(SS(i):EE(i)));
            if i~=length(SS)
                hover1(i,2)=sum(DiffLoc(EE(i):SS(i+1)-1));
                hover1(i,3)=range(Loc(EE(i):SS(i+1)));
            else
                if EE(end)>length(Mark)
                    hover1(i,:)=[];
                else
                    hover1(i,2)=sum(DiffLoc(EE(i):length(Mark)));
                    hover1(i,3)=range(Loc(EE(i):length(Mark)+1));
                end
            end
        end  
        if SS(1)>1
            hover1=[[SS(1)-1 sum(DiffLoc(1:SS(1)-1)) range(Loc(1:SS(1))) SeqNo];hover1];                
        end
    end
    CandStepSize(SeqNo,:)=[size(hover1,1) size(move1,1) range(Loc)/length(Mark) mean(DiffLoc)];    
    % CandStepSize: hover step number; jump step number; predict step size
   % based on range (earlier version); predict step size based on sum
    hover = [hover; hover1];
    move = [move; move1];
    
    % duration (time bins), total travel pos bins(sum), total travel pos bins(range); SeqNo
end
% disp('wait')

eval([label  '_hover = hover;'])
eval([label  '_move = move;'])
eval([label  '_CandStepSize = CandStepSize;'])
eval([label  '_seqisjump_bin = seqisjump_bin;'])
eval([label  '_seqisjump_btwn = seqisjump_btwn;'])
eval([label  '_seqstepst_btwn = seqstepst_btwn;'])
eval([label  '_seqtimes_btwn = seqtimes_btwn;'])
eval([label  '_seqtimes_bin = seqtimes_bin;'])

eval([label  '_seqindex_btwn = seqindex_btwn;'])
eval([label  '_seqindex_bin = seqindex_bin;'])

 
save(thisdir, [label  '_hover'],...
    [label  '_move'],[label  '_CandStepSize'],...
    [label  '_seqisjump_bin'],[label '_seqisjump_btwn'],[label  '_seqstepst_btwn'],...
    [label  '_seqtimes_btwn'],[label  '_seqtimes_bin'],...
    [label  '_seqindex_btwn'],[label  '_seqindex_bin'],'-append')