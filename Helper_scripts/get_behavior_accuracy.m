function armacc = get_behavior_accuracy(thisdir)

load(thisdir,'laps_singlepass','armpos','rundat')
[~,rn] = max(rundat,[],2);
armacc = [];
for irun = 1:max(rn)
    il = min(laps_singlepass(rn==irun));
    currarm = armpos(find(laps_singlepass==il,1,'first')); % 3 is center
    nextarm = armpos(find(laps_singlepass==il,1,'last'));
    wenttocenter = true;
    wenttoarm = false;

    if currarm~=3

        wenttoarm = true;
        if nextarm==3
            lastarm = currarm;            
            armacc = [armacc;currarm nextarm true];
        else
            lastarm = nextarm;            
            armacc = [armacc;currarm nextarm false];
        end
    elseif currarm==3
        armacc = [armacc;currarm nextarm NaN];
        lastarm = nextarm;
    end
    

    for ilap = il+1:max(laps_singlepass(rn==irun))
        currarm = armpos(find(laps_singlepass==ilap,1,'first')); % 3 is center
        nextarm = armpos(find(laps_singlepass==ilap,1,'last'));
        if nextarm==3             
           wenttocenter = true;
           wenttoarm = false;
           armacc = [armacc;currarm nextarm true];
        elseif nextarm==1 || nextarm==2
            if nextarm~=lastarm && wenttocenter && ~wenttoarm
                armacc = [armacc;currarm nextarm true];
                lastarm = nextarm;
            else
                armacc = [armacc;currarm nextarm false];
            end
            wenttoarm=true;
            wenttocenter = false;                   
        end
            

    end
    
end
armacc(:,4) = armacc(:,3);
% armacc(armacc(:,2)==3,4) = NaN; %ABL changed 10/12/2020 from this
armacc(armacc(:,1)~=3,4) = NaN; %ABL changed 10/12/2020 to this
%arm started on, arm went to, rewarded, alternation accuracy