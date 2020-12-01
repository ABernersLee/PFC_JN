
allcells1 = [];
allcells2 = [];
smi = [];
cd('D:\XY_matdata\AllSessions')
d2 = dir('*.mat');
am1 = []; am2 = []; amC1 = []; amC2 = [];
for id =  1:size(d2,1)
    thisdir = d2(id).name;
   
    [daycells1,am] =  get_PFC_hoverjump_slopeinout(thisdir,'RP');
    allcells1 = cat(1,allcells1,daycells1);
    am1 = cat(1,am1,am);
    
   
    [daycells2,am] =  get_PFC_hoverjump_slopeinout(thisdir,'SD');
    allcells2 = cat(1,allcells2,daycells2);
    am2 = cat(1,am2,am);
            
    load(thisdir,'RP_Cand_sig_modu_include')
    if ~exist('RP_Cand_sig_modu_include','var')
        sig_modu_include = make_PFC_candeventtriggeredmat(thisdir,'RP');
    else
        sig_modu_include = RP_Cand_sig_modu_include;
        clear RP_Cand_sig_modu_include
    end
    smi = cat(1,smi,sig_modu_include);
    
    disp(['done ' num2str(id)])
end

r1 = corr([1:10]',allcells1(smi(:,2)>0 & smi(:,1)==1 & smi(:,3)==1,:,1)');
r2 = corr([1:10]',allcells1(smi(:,2)>0 & smi(:,1)==1 & smi(:,3)==1,:,2)');
r3 = corr([1:10]',allcells1(smi(:,2)<0 & smi(:,1)==1 & smi(:,3)==1,:,1)');
r4 = corr([1:10]',allcells1(smi(:,2)<0 & smi(:,1)==1 & smi(:,3)==1,:,2)');
ranksum(r1-r2,0)
ranksum(r3-r4,0)
ranksum(r1-r2,r3-r4)

if 1 
mmlab = {'All Cells';'Neg Mod';'Pos Mod';'Sig Neg Mod';'Sig Pos Mod'};
rplab = {'Ripple Cand Events';'SD Cand Events'};
for irpsd = 1:2
    if irpsd == 1
        dat1 = allcells1;
    elseif irpsd == 2
        dat1 = allcells2;
    end
    for mm = 2:5    
        if mm == 1
            ind = smi(:,3)==1;
        elseif mm == 2 %neg mod
            ind = smi(:,3)==1 & smi(:,2)<0;
        elseif mm == 3 %pos mod
            ind = smi(:,3)==1 & smi(:,2)>0;
        elseif mm == 4 % sig neg mod
            ind = smi(:,3)==1 & smi(:,2)<0 & smi(:,1)==1;
        elseif mm==5 %sig pos mod
            ind = smi(:,3)==1 & smi(:,2)>0 & smi(:,1)==1;
        end
        ind2 = 1:size(dat1,2);
        
        
        
        figure; hold on
        
        
        for ii = 1:2
            dat = dat1(:,:,ii);
            r1 = corr([1:size(ind2,2)]',dat(ind,ind2)');
            subplot(1,2,ii)
            sem = nanstd(zscore(dat(ind,ind2)'),[],2)./sqrt(size(dat(ind,ind2),1));
            meandat = nanmean(zscore(dat(ind,ind2)'),2);
            rev = meandat-sem;
            plot(meandat,'k','LineWidth',2)
            patch([1:size(ind2,2) size(ind2,2):-1:1],[meandat+sem;rev(end:-1:1)],'black','FaceAlpha',.3,'EdgeAlpha',0)
            [r2,p2] = corr([1:size(ind2,2)]',meandat);
            tt = title([mmlab{mm} ' ' rplab{irpsd} ' r = ' num2str(r2) ' p = ' num2str(p2)]);
            if p2<.05
                tt.Color = 'r';
            end
        end
    end
end
disp('wait')
end