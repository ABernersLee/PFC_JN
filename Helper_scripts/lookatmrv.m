
allcells1 = [];
allcells2 = [];
allcellsA = [];
allcellsB = [];
smi = [];
cd('D:\XY_matdata\AllSessions')
% std_cutoff = 2;
% transferdata
% transferLFPandTT
d2 = dir('*.mat');
spikecutoff = 4;
for id =  1:size(d2,1)
    thisdir = d2(id).name;
   
    [daycells1,daycellsA] = get_PFC_hoverjump_spikeprob_mrv(thisdir,'RP',spikecutoff);
    if ~isempty(daycells1)
        allcells1 = cat(1,allcells1,daycells1);
    end
    allcellsA = cat(1,allcellsA,daycellsA);
   
%     [daycells2,daycellsB] = get_PFC_hoverjump_spikeprob_mrv(thisdir,'SD',spikecutoff);
%     allcells2 = cat(1,allcells2,daycells2);
%     allcellsB = cat(1,allcellsB,daycellsB);
        
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

%%
mmlab = {'All Cells';'Neg Mod';'Pos Mod';'Sig Neg Mod';'Sig Pos Mod'};
rplab = {'Ripple Cand Events';'SD Cand Events'};
for irpsd = 1:2
    if irpsd == 1
        dat = allcells1;
    elseif irpsd == 2
        dat = allcells2;
    end
    for mm = 4:5    
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
        ind2 = 1:size(dat,2);
        
        r1 = corr([1:size(ind2,2)]',dat(ind,ind2)');
        nanmean(r1);
        figure; hold on
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
disp('wait')

%%

