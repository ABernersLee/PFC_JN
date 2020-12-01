
%quick look at gamma modulation
% 
% load(thisdir,[label  '_seqisjumpbin'],[label  '_seqstepst'])
% eval(['isjump = ' label  '_seqisjumpbin;'])
% eval(['isjumpst = ' label  '_seqstepst;'])
% 
% sm = 15;
% [~,~,i] = histcounts(seqphases,0:(360/numbins):360);
% clear meanss; for ideg = 1:numbins; meanss(ideg,1) = mean(seqHPspikes(i==ideg)); end
% output = smoothts([meanss;meanss;meanss;meanss]','b',sm);
% o1 = output(length(meanss)+1:(3*length(meanss)));
% 
% clear meanss; for ideg = 1:numbins; meanss(ideg,1) = mean(isjump(i==ideg)==0); end
% output = smoothts([meanss;meanss;meanss;meanss]','b',sm);
% o2 = output(length(meanss)+1:(3*length(meanss)));
% 
% clear meanss; for ideg = 1:numbins; meanss(ideg,1) = mean(seqPFCspikes(i==ideg)); end
% output = smoothts([meanss;meanss;meanss;meanss]','b',sm);
% o3 = output(length(meanss)+1:(3*length(meanss)));
% 
% clear meanss; for ideg = 1:numbins; meanss(ideg,1) = mean(isjumpst(i==ideg)==1); end
% output = smoothts([meanss;meanss;meanss;meanss]','b',sm);
% o4 = output(length(meanss)+1:(3*length(meanss)));
% 
% clear meanss; for ideg = 1:numbins; meanss(ideg,1) = mean(isjumpst(i==ideg)==-1); end
% output = smoothts([meanss;meanss;meanss;meanss]','b',sm);
% o5 = output(length(meanss)+1:(3*length(meanss)));
% 
% figure; hold on, 
% plot(zscore(o1),'k','LineWidth',3), 
% plot(zscore(o2),'r','LineWidth',3); 
% plot(zscore(o3),'g','LineWidth',3);
% plot(zscore(o4),'b','LineWidth',3);
% % plot(zscore(o5),'c','LineWidth',3);