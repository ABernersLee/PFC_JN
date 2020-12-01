%exploring strategies of the rats in early learning

pe = perms(3:-1:1); pe = pe(:,1:2);
hs = NaN(size(d2,1),1);
s = NaN(size(d2,1),4);
half = NaN(size(d2,1),3); halfAlt = half;
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    armacc = get_behavior_accuracy(thisdir);   
    for ip = 1:size(pe,1)
        armacc(armacc(:,1)==pe(ip,1) & armacc(:,2)==pe(ip,2),5) = ip;
    end
    h = histc(armacc(:,5),1:size(pe,1));
    rl = zeros(3,3);
    
    for ip = 1:size(pe,1)
        rl(pe(ip,2),pe(ip,1)) = h(ip)./sum(h);
    end
    %arm started on, arm went to, rewarded, alternation accuracy
%     figure; imagesc(rl); axis xy; colormap gray
%     hold on; xlabel('From'), ylabel('To')
%     title([num2str(id)])
    h = histc([armacc(:,1);armacc(end,2)],1:3);
    hs(id,:) = min(h./sum(h))*100;
    
    recent = NaN(size(armacc,1),1);
    total = recent;
    notrecent = recent; last = recent;
    for itrial = 2:size(armacc,1)
        notrecent(itrial) = setdiff(1:3,armacc(itrial-1:itrial,1));
        last(itrial) = armacc(itrial-1,1);
        if sum(armacc(1:itrial-2,3)==1)>0 && itrial>2
            if sum(armacc(1:itrial-2,3)==1 & armacc(1:itrial-2,2)~=armacc(itrial,1))>0
            recent(itrial) = armacc(find(armacc(1:itrial-2,3)==1 & armacc(1:itrial-2,2)~=armacc(itrial,1),1,'last'),2);
            end
            total1 = NaN(3,1);
            for iarm = 1:3
                total1(iarm) = sum(armacc(armacc(:,2)==iarm,3))./sum(armacc(:,2)==iarm);
            end
            total1(armacc(itrial,1)) = 0;
            if max(total1)>0
                [~,total(itrial)] = max(total1);
            end
        end
    end
    t = [recent total notrecent last];
    t(sum(isnan(t),2)>0,:) = NaN; %checking to see if you restrict it to trials where there is a difference between the arms in all measures, recently rewarded still wins out
    s(id,:) = (sum(t==armacc(:,2))./sum(~isnan(t)))*100;
    half(id,:) = [nanmean(armacc(1:floor(size(armacc,1)/2),3)) nanmean(armacc(end-(floor(size(armacc,1)/2)-1):end,3)) nanmean(armacc(:,3))];
    halfAlt(id,:) = [nanmean(armacc(1:floor(size(armacc,1)/2),4)) nanmean(armacc(end-(floor(size(armacc,1)/2)-1):end,4)) nanmean(armacc(:,4))];
    
    
%     half(id,:) = [nanmean(armacc(1:3,3)) nanmean(armacc(4:end,3)) nanmean(armacc(:,3))];
%     halfAlt(id,:) = [nanmean(armacc(1:3,4)) nanmean(armacc(4:end,4)) nanmean(armacc(:,4))];
end
% mean(s)
% median(s)
%as a whole, rats are more likely to visit the most recently rewarded arm
%than the one with the highest current reward rate, the arm not visited
%recently or the arm visited more recently

%lots of variability
[mm,m] = max(s,[],2);
m(sum(mm==s,2)>1) = NaN;
histc(m,1:4);

mx= mm==s;
sum(mx);
mx = (1./sum(mx,2))./mx;
mx(isinf(mx)) = 0;
sum(mx);
%most frequent type of session is to have the rat choose the arm with the
%highest total reward rate
sum(s>50)
% sum(s(:,1)>50 | s(:,2)>50)
%in most of the sessions (7/11), the rat was more likely to choose the arm with
%the highest current reward rate, given its past choices. In those sessions,
%rats chose the arm with the highest reward rate 66% of the time
%(chance is 50%), and across all sessions rats chose it
%54% of the time.
%In 3 of the other 4 sessions, rats were more likely to
%choose an arm if it had not been visited recently (61% of the time), and
%in one session the rat tended to visit the arm he had most recently
%visited (75% of the time) but also was performing the alternation rule 56% of the time.
sum(s(s(:,2)<=50,:)>50)
mean(s(s(:,2)>50,2))
mean(s(s(:,2)<=50 & s(:,3)>50,3))

% look at across days whether the reward rates and alternation rates
% improved
across = [1 3; 4 7; 8 10];
mean(halfAlt(across(:,1),3))
mean(halfAlt(across(:,2),3))
mean(half)
mean(halfAlt)
days = [1:3 1:4 1:3 NaN];
[r,p] = corr(half(:,3),days','rows','complete','type','Kendall');
[r,p] = corr(halfAlt(:,3),days','rows','complete','type','Kendall');