cd('D:\XY_matdata\AllSessions')
aa1 = []; aa2 = []; smi = []; aa3 = []; tt1 = [];
d2 = dir('*.mat');
corrcut = 0; jumpcut = 1.1;
label = 'RP';
for id =  1:size(d2,1)
    thisdir = d2(id).name;
%     load(thisdir,[label '_pSSDarm'],[label '_moduarm'])
%     if ~exist([label '_pSSDarm'],'var')
%         get_PFC_armtriggered_modusig(thisdir,label)
%     end
    [armmodu,armmodudur,tot] = make_durationmodu(thisdir,label,corrcut,jumpcut);
    aa1 = cat(1,aa1,armmodu);
    aa2 = cat(1,aa2,armmodudur);
    tt1 = cat(1,tt1,tot);
    
    rss = make_duration_dur(thisdir,label,corrcut,jumpcut);
    aa3 = cat(1,aa3,rss);
    
    load(thisdir,'RP_Cand_sig_modu_include')
    if ~exist('RP_Cand_sig_modu_include','var')
        sig_modu_include = make_PFC_candeventtriggeredmat(thisdir,'RP');
    else
        sig_modu_include = RP_Cand_sig_modu_include;
        clear RP_Cand_sig_modu_include
    end
    smi = cat(1,smi,sig_modu_include);    
end
%%
ind = smi(:,2)>0 & smi(:,1)==1 & smi(:,3)==1;
dat1 = aa1(ind,:);
dat2 = aa2(ind,:);
[r1,p1] = corr(dat1(:),dat2(:),'rows','complete')

ind = smi(:,2)<0 & smi(:,1)==1 & smi(:,3)==1;
dat1 = aa1(ind,:);
dat2 = aa2(ind,:);
[r2,p2] = corr(dat1(:),dat2(:),'rows','complete')

ind = smi(:,1)==1 & smi(:,3)==1;
dat1 = aa1(ind,:);
dat2 = aa2(ind,:);
[r3,p3] = corr(dat1(:),dat2(:),'rows','complete')


