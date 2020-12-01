function prospective_coding(thisdir)
% warning ('off','all')
% warning

load(thisdir,'other_cells','spikedata')
load(thisdir, 'armpos', 'behave_change_log', 'headingarm','pos','dirname','laps_singlepass','vel','dirdat')

nS = 1000;

velcutoff = 10;
L = laps_singlepass;
Labels = cell(3,1);
sd = spikedata(ismember(spikedata(:,2),other_cells),:);
FRmap = NaN(length(other_cells),max(L));
for ilap = 1:max(L)        
    currarm = armpos(find(L==ilap,1,'last'));
    lastarm = armpos(find(L==ilap & armpos==currarm,1,'first')-1);        
    ind = armpos==lastarm & L == ilap & headingarm==currarm & vel>velcutoff & dirdat==0;  
    if sum(ind)>0
        Labels{lastarm} = cat(1,Labels{lastarm},[ilap currarm]);
        spkind = ismember(sd(:,3),find(ind));
        [f,~,~] = histcounts(sd(spkind,2),[other_cells-.5;other_cells(end)]);
        t = sum(ind)/30;
        FRmap(:,ilap) = f./t;
    end
end

acc = []; accS = [];
for iarm = 1:3
    Label = Labels{iarm}(:,2);
    FR = FRmap(:,Labels{iarm}(:,1));
    [Train1, ~] = crossvalind('LeaveMOut', length(Label), 1);
    if length(Label)>3
    for j = 1:length(Label)           
        train = circshift(Train1,1);
        test = ~train;           
%             guess = classify(FR(:,test)',FR(:,train)',Label(train),'linear');        
        LL = Label(train);
        
        Mdl = fitcdiscr(FR(:,train)',Label(train),'discrimType','pseudoLinear');
        guess = predict(Mdl,FR(:,test)');
        %         CoeffsPFC1(PFCind,:) = Mdl.Coeffs(1,2).Linear;
        
%         B = mnrfit(FR(:,train)',LL);
%         f = mnrval(B,FR(:,test)');
%         [~,guess] = max(f,[],2);

        acc = cat(1,acc,guess==Label(test));
        gs = NaN(nS,1);
        for ishuff = 1:nS
%             B = mnrfit(FR(:,train)',LL(randperm(length(LL))));
%             f = mnrval(B,FR(:,test)');
%             [~,sguess] = max(f,[],2);            
            
              Mdl = fitcdiscr(FR(:,train)',LL(randperm(length(LL))),'discrimType','pseudoLinear');
%                 CoeffsPFC1(PFCind,i) = Mdl.Coeffs(1,2).Linear;
              sguess = predict(Mdl,FR(:,test)');
              gs(ishuff,1) = sguess==Label(test);
        end
        accS = cat(2,accS,gs);
    end
    end
end
p = (sum(mean(accS,2)>=mean(acc))+1)/(size(accS,1)+1);
% prospectivecoding_p = ranksum(acc,accS,'tail','right');
% save(thisdir,'prospectivecoding_p','-append')
% warning ('on','all')
% warning('query','all')
