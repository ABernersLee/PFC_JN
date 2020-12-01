function make_singleandjointreplay_fromPosterior(thisdir,label,numshuff)
disp(['Start make_singleandjointreplay_fromPosterior ' label])
load(thisdir,[label '_InMatrix'],[label '_OutMatrix'],[label '_Index'],'armpos',[label '_CandEventTimes'],'armposindex')
% load('Ripple_Events')
% Cand = [Ripple_Events(:,1)-.02 Ripple_Events(:,2)+.02];

[~,armpos] = max(armposindex,[],2);
eval(['InMatrix = ' label  '_InMatrix;'])
eval(['OutMatrix = ' label  '_OutMatrix;'])
eval(['Index = ' label  '_Index;'])
eval(['Cand = ' label  '_CandEventTimes(:,1:2);'])
Cand = [Cand(:,1)-.02 Cand(:,2)+.02];

EstBin = .02;
N=round(diff(Cand')/EstBin)';
t=find(mod(N,2));
Cand(t,:)=[Cand(t,1)/2+Cand(t,2)/2-EstBin*floor(N(t)/2)-EstBin/2 Cand(t,1)/2+Cand(t,2)/2+EstBin*floor(N(t)/2)+EstBin/2];
t=find(~mod(N,2));
Cand(t,:)=[Cand(t,1)/2+Cand(t,2)/2-EstBin*N(t)/2 Cand(t,1)/2+Cand(t,2)/2+EstBin*N(t)/2];
clear N t EstBin


D=OutMatrix+InMatrix;   
% if ilab == 1
%     WCcutoff = .5; %was .3, changed 2/4/19 to check, 2/6/19 changed back to make more flexible
% elseif ilab == 2
WCcutoff = .3; 
% end
timecutoff = round(.05/.005); %to make it 50ms, changed to 40ms, changed back 2/4/19
MAPcutoff = 4; %new @ 4, was 5, was 4, changed 2/4/19, 2/6/19 changed back
armcovcutoff = .5;
JDcutoff = 1; % was .7, changed 2/4/19 to check
gauss_x=-20:1:20;
gauss_y=normpdf(gauss_x,0,1);
gaussFilter=gauss_y/sum(gauss_y);

replayqual = [];
replayqual.currentnumber = 0;
replayqual.singlebothjoint = []; replayqual.times = []; replayqual.post= []; replayqual.postseqindex = [];
replayqual.wc = [];replayqual.armcov = []; replayqual.distance = []; replayqual.maxjump = []; replayqual.dirbias = [];
replayqual.replayarm = [];replayqual.jointreplayarm = [];replayqual.jointjumptime = [];
replayqual.seq = []; replayqual.seqtimes = []; replayqual.seqreplayindex = [];
replayqual.seqjointjump = []; replayqual.stnd = []; replayqual.seqmax = []; replayqual.shuffle_p = [];

CandEventAssign = NaN(length(Index),1); %should have one fewer row, but i fix it later instead
tic
for c = 1:length(Index)-1    
    replaybins = 4*Index(c)+1:4*Index(c+1)-3;
    times = Cand(c,1):.005:Cand(c,2)-.005;
    MAP_arms = NaN(length(replaybins),max(armpos));                 
    
    for iarm = 1:max(armpos)    
        segment = D(replaybins,armpos==iarm); 
        MAP1 = max(segment,[],2);          
        withwings = conv(MAP1, gaussFilter);
        MAP=withwings(((size(gauss_x,2)-1)/2)+1:end-(size(gauss_x,2)-1)/2);
        MAP_arms(:,iarm) = MAP; 
    end    
    cutoff = ones(size(MAP_arms)).*(1/length(armpos))*MAPcutoff; 
    cut = MAP_arms>cutoff;
    cutdiff = [zeros(1,max(armpos));diff(cut)];
    ireplays = []; iis = [];
    
    %evaluate if any arms pass criterion of having a high posterior for at least 50ms
    for iarm = 1:max(armpos)
        if sum(cut(:,iarm))>=timecutoff
            if cut(1,iarm)==1
                stt = [1;find(cutdiff(:,iarm)==1)];
            else
            stt = find(cutdiff(:,iarm)==1);
            end
            endd = find(cutdiff(:,iarm)==-1)-1;
            if size(endd,1)>size(stt,1)
                iis(:,1) = [1;stt];
            else iis(:,1) = stt;
            end
            if size(endd,1)<size(stt,1)
                iis(:,2) = [endd;size(cutdiff,1)];                
            else iis(:,2) = endd;
            end                            
            iis((iis(:,2)-iis(:,1))<timecutoff,:) = [];
            
            if ~isempty(iis)
                sb = find([false;(iis(2:end,1)-iis(1:end-1,2))<timecutoff+1]);
                if ~isempty(sb)
                    for i = 1:length(sb)
                        iis(sb(i)-1,2) = iis(sb(i),2);
                        iis(sb(i),:) = [NaN NaN];
                    end
                    iis(isnan(iis(:,1)),:) = []; 
                end
                
                armcovs = NaN(size(iis,1),1); wcs = armcovs; mjump = armcovs;
                for i = 1:size(iis,1)
                    Mat = D(replaybins(iis(i,1):iis(i,2)),armpos==iarm)';
                    I2=sum([1:size(Mat,1)]'*ones(1,size(Mat,2)).*Mat);
                    
                    segment = Mat';
                    extend_x = repmat(1:size(segment,1),size(segment,2),1);
                    extend_y = repmat([1:size(segment,2)]',1,size(segment,1));
                    wcs(i) = abs(weighted_corr(extend_x,extend_y,segment'));   
%                     [x,y] = find(segment==repmat(max(segment,[],2),[1 size(segment,2)]));
%                     seqmax = NaN(size(I2));
%                     seqmax(x) = y;         
                    
                    seqmax = [1:size(Mat,1)]*(Mat./sum(Mat,1));
                    
                    armcovs(i) = range(seqmax)/sum(armpos==iarm);
                    mjump(i) = max(abs(diff(seqmax)))/sum(armpos==iarm);
                end
                
                ireplays = [ireplays; iis iarm*ones(size(iis,1),1) wcs armcovs mjump];                                
            end
        end       
        iis = [];
    end
    
    
    allpassedtime = ireplays;
    if ~isempty(ireplays)
%         ireplays = ireplays(ireplays(:,4)>WCcutoff  & ireplays(:,5)>armcovcutoff & [(abs(diff(ireplays(:,1:2)'))*.005)>=.049]',:);
        ireplays = ireplays(ireplays(:,4)>WCcutoff & ireplays(:,6)<JDcutoff & ireplays(:,5)>armcovcutoff & [(abs(diff(ireplays(:,1:2)'))*.005)>=.05]',:);
    end
    joi = [];
    
    %evaluate any possible joint replays
    if size(ireplays,1)>1
        if length(unique(ireplays(:,3)))>1              
             jj = combnk(1:size(ireplays,1),2);   
             jarm = combnk(ireplays(:,3),2);  
             
             jwcs = zeros(size(jj,1),1); jjd = jwcs;
             jind_all = NaN(size(jj,1),2);  
             jjind = NaN(size(jj,1),2,2);
             jord_all = cell(size(jj,1),1);  
             for jr = 1:size(jj,1)    
                 if diff(jarm(jr,:))==0
                     continue
                 end
                    % ABL added 10/11/18
                   [~,frstarm] = min([ireplays(jj(jr,1),1) ireplays(jj(jr,2),1)]);
                   if frstarm==2
                       jarmsave = jarm(jr,1);
                       jarm(jr,1) = jarm(jr,2);
                       jarm(jr,2) = jarmsave;
                       jjsave = jj(jr,1);
                       jj(jr,1) = jj(jr,2);
                       jj(jr,2) = jjsave;
                   end
                   
                   fwd1 = find(armpos==ireplays(jj(jr,1),3));% to flip the arm around
                   fwd2 = find(armpos==ireplays(jj(jr,2),3)); %the first one needs to flip around going from high to low, the second one needs to continue, going low to high                   
                   jord = [fwd1(end:-1:1);fwd2]; % to be facing eachother (middles in middle)                                      
                   jsegord = D(replaybins,jord);
                   jind = [min(ireplays(jj(jr,1),1)) max(ireplays(jj(jr,2),2))];
                   jjind(jr,:,:) = [ireplays(jj(jr,1),1) ireplays(jj(jr,2),1);ireplays(jj(jr,1),2) ireplays(jj(jr,2),2)];
                   jseg = jsegord(jind(1):jind(2),:); %using time bins from the beginning of the first to the end of the last
                   extend_x = repmat(1:size(jseg,1),size(jseg,2),1);
                   extend_y = repmat([1:size(jseg,2)]',1,size(jseg,1));
                   jwc = weighted_corr(extend_x,extend_y,jseg'); 
                   [x,y] = find(jseg'==repmat(max(jseg',[],1),[size(jseg',1) 1]));
                   seqmax = NaN(size(I2));
                   seqmax(y) = x;                                       
                   jwcs(jr,1) = abs(jwc);            
                   jjd(jr,1) = max(abs(diff(seqmax)))./size(jseg,2);
                   jind_all(jr,:) = jind;
                   jord_all{jr} = jord;
             end          
             joi = find(abs(jwcs)>=WCcutoff); %ABL changed 10/11/18 from % & jjd<=JDcutoff);
             if length(joi)>1
                   joi = joi(jwcs(joi)==max(jwcs(joi)));
                   joi = joi(1);
             end
        end
    end
    
    %Now add any single or joint replays to the matricies of replay qualities
    if ~isempty(joi)
        %add as a joint replay and also as two single replays                 
        CandEventAssign(c) = 2;
        j = mean([min(jjind(joi,:,2)) max(jjind(joi,:,1))]);
        jj = [mean(diff(times))*(j-1)+times(1) round(j)];        
        replayqual = addreplay(D,jind_all(joi,1),jind_all(joi,2),jord_all{joi},3,replayqual,jarm(joi,:),InMatrix,OutMatrix,times,jj,replaybins,numshuff);
        for ijarm = 1:2            
            thisarm = jarm(joi,ijarm);         
            % add both singles
            replayqual = addreplay(D,jjind(joi,1,ijarm),...
                jjind(joi,2,ijarm),find(armpos==thisarm),2,replayqual,jarm(joi,ijarm),InMatrix,OutMatrix,times,NaN,replaybins,numshuff);         
        end
    else
        if ~isempty(ireplays)
            CandEventAssign(c) = 1;
            if size(ireplays,1)==1
                %add the only single replay
                replayqual = addreplay(D,ireplays(1),ireplays(2),find(armpos==ireplays(3)),1,replayqual,ireplays(3),InMatrix,OutMatrix,times,NaN,replaybins,numshuff);               
            elseif size(ireplays,1)>1
                %add best of the single replays
                bestarm = find(max(abs(ireplays(:,4)))==abs(ireplays(:,4)));
                replayqual = addreplay(D,ireplays(bestarm,1),ireplays(bestarm,2),find(armpos==ireplays(bestarm,3)),1,replayqual,ireplays(bestarm,3),InMatrix,OutMatrix,times,NaN,replaybins,numshuff);                
            end
        elseif ~isempty(allpassedtime)            
            %add best of the bad replays
             if size(allpassedtime,1)==1
                %add the only bad replay
                CandEventAssign(c) = 1;
                replayqual = addreplay(D,allpassedtime(1),allpassedtime(2),find(armpos==allpassedtime(3)),0,replayqual,allpassedtime(3),InMatrix,OutMatrix,times,NaN,replaybins,numshuff);           
            elseif size(allpassedtime,1)>1
%                 add best of the bad replays if one is best, changed 2/5/19       
                bestarm = find(max(abs(allpassedtime(:,4)))==abs(allpassedtime(:,4)));
                if length(bestarm)==1
                    CandEventAssign(c) = 1;
                    replayqual = addreplay(D,allpassedtime(bestarm,1),allpassedtime(bestarm,2),find(armpos==allpassedtime(bestarm,3)),0,replayqual,allpassedtime(bestarm,3),InMatrix,OutMatrix,times,NaN,replaybins,numshuff);             
                else
                   CandEventAssign(c) = 0; 
                end
            end
        else
            % get rid of this candidate event
            CandEventAssign(c) = 0;
        end
        
    end   
    if rem(c,500)==0
        t = toc;
       disp([num2str(c) ' Out of ' num2str(length(Index)-1) ' Events Finished in ' num2str(t/60) ' Minutes']) 
       tic
    end
end

%break down replayqual into seperate matriceis
% disp('wait')
replayqual.postseqindex = [1; replayqual.postseqindex(1:end-1)+1];

eval([label  '_replay_singlebothjoint =replayqual.singlebothjoint;'])
eval([label  '_replay_times =replayqual.times;'])
eval([label  '_replay_jointjumptime =replayqual.jointjumptime;'])
eval([label  '_replay_post =replayqual.post;'])
eval([label  '_replay_postseqindex =replayqual.postseqindex;'])
eval([label  '_replay_corr =replayqual.wc;'])
eval([label  '_replay_armcov =replayqual.armcov;'])
eval([label  '_replay_distance =replayqual.distance;'])
eval([label  '_replay_maxjump =replayqual.maxjump;'])
eval([label  '_replay_dirbias =replayqual.dirbias;'])
eval([label  '_replay_replayarm =replayqual.replayarm;'])
eval([label  '_replay_jointreplayarm =replayqual.jointreplayarm;'])
eval([label  '_replay_seqmean =replayqual.seq;'])
eval([label  '_replay_seqmax =replayqual.seqmax;'])
eval([label  '_replay_seqtimes =replayqual.seqtimes;'])
eval([label  '_replay_seqreplayindex =replayqual.seqreplayindex;'])
eval([label '_CandEventAssign = CandEventAssign;'])
eval([label  '_replay_seqjointjump =replayqual.seqjointjump;'])
eval([label  '_replay_stnd =replayqual.stnd;'])
eval([label  '_replay_shuffle_p =replayqual.shuffle_p;'])


% look = replayqual.jointreplayarm(replayqual.singlebothjoint==3 & replayqual.shuffle_p<.05,:);
% confusionmat(look(:,1),look(:,2))
[sum(replayqual.singlebothjoint==3 & replayqual.shuffle_p<.05 & abs(replayqual.wc)>WCcutoff) sum(replayqual.singlebothjoint==1 & replayqual.shuffle_p<.05 & abs(replayqual.wc)>WCcutoff & replayqual.maxjump<JDcutoff)]


save(thisdir,[label  '_replay_singlebothjoint'],[label '_CandEventAssign'],[label  '_replay_times'],...
    [label  '_replay_jointjumptime'],[label  '_replay_seqmean'],...
    [label  '_replay_post'],[label  '_replay_postseqindex'],[label  '_replay_seqjointjump'],...
    [label  '_replay_corr'],[label  '_replay_armcov'],[label  '_replay_distance'],[label  '_replay_maxjump'], ...
    [label  '_replay_dirbias'],[label  '_replay_replayarm'],[label  '_replay_jointreplayarm'],[label  '_replay_seqmax'], ...
    [label  '_replay_seqtimes'],[label  '_replay_seqreplayindex'],[label  '_replay_stnd'],[label  '_replay_shuffle_p'],'-append')

disp(['Done with make_singleandjointreplay_fromPosterior ' label])

