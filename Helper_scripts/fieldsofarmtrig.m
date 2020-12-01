cd('D:\XY_matdata\AllSessions')
d2 = dir('*.mat');

ilabel = 1;
if ilabel == 1
    label = 'RP';
elseif ilabel == 2
    label = 'SD';
end
M = []; pSD = []; smi = []; AllIn = []; AllOut = [];
dirbias = [];
for idir = 1:20
    thisdir = d2(idir).name;    
    load(thisdir,[label '_pSSDarm'],[label '_moduarm'])
    if ~exist([label '_moduarm'],'var')
        [p_SSD,modu] = get_PFC_armtriggered_modusig(thisdir,label);
    else
        eval(['p_SSD = ' label '_pSSDarm;'])
        eval(['modu = ' label '_moduarm;'])
        clear([label '_pSSDarm'],[label '_moduarm'])
    end
    pSD = cat(1,pSD,p_SSD);
    disp(['Day ' num2str(idir) ' ' num2str((sum(p_SSD<.05)./sum(~isnan(p_SSD)))*100)])
    M = cat(1,M,modu');
    load(thisdir,'OpenInFR','OpenOutFR')
    if ~exist('OpenInFR','var')
        [OpenInFR,OpenOutFR] = make_directionalfields_open(thisdir);
    end
    load(thisdir,'other_cells')
%     load(thisdir,'OutFR','InFR','other_cells')
%     InPFC = InFR(other_cells,:); clear InFR
%     OutPFC = OutFR(other_cells,:); clear OutFR 
%     dirbias = cat(1,dirbias,(sum(InPFC,2)-sum(OutPFC,2))./(sum(InPFC,2)+sum(OutPFC,2)));    
%     disp(['Done with ' num2str(idir)])
    
    load(thisdir,[label '_Cand_sig_modu_include'])
    if ~exist([label '_Cand_sig_modu_include'],'var')
        sig_modu_include = make_PFC_candeventtriggeredmat(thisdir,'RP');
    else
        sig_modu_include = RP_Cand_sig_modu_include;
        clear RP_Cand_sig_modu_include
    end
    smi = cat(1,smi,sig_modu_include);
%     if size(OpenOutFR,1)<size(AllOut,1)
%         OpenOutFR = cat(1,OpenOutFR,NaN(size(AllOut,1)-size(OpenOutFR,1),size(OpenOutFR,2),size(OpenOutFR,3)));
%         OpenInFR = cat(1,OpenInFR,NaN(size(AllOut,1)-size(OpenInFR,1),size(OpenInFR,2),size(OpenInFR,3)));
%     elseif size(OpenOutFR,2)>size(AllOut,2) && idir>1
%         AllOut = cat(2,AllOut,NaN(size(AllOut,1),size(OpenOutFR,2)-size(AllOut,2),size(AllOut,3)));
%         AllIn = cat(2,AllIn,NaN(size(AllIn,1),size(OpenInFR,2)-size(AllIn,2),size(AllIn,3)));
%     end
%     AllOut = cat(3,AllOut,OpenOutFR(:,:,other_cells));
%     AllIn = cat(3,AllIn,OpenInFR(:,:,other_cells));    
    AllOut = OpenOutFR(:,:,other_cells);
    AllIn = OpenInFR(:,:,other_cells);  
    AllAll = AllIn+AllOut;
    InAll2 = AllIn./repmat(max(max(AllIn,[],2),[],1),[size(AllIn,1) size(AllIn,2) 1]);
    OutAll2 = AllOut./repmat(max(max(AllOut,[],2),[],1),[size(AllOut,1) size(AllOut,2) 1]);
    AllAll2 = AllAll./repmat(max(max(AllAll,[],2),[],1),[size(AllAll,1) size(AllAll,2) 1]);
    clear OpenInFR OpenOutFR
    figure; hold on,subplot(1,3,1); imagesc(mean(InAll2(:,:,p_SSD<.05),3)-mean(InAll2(:,:,p_SSD>.05),3))
    subplot(1,3,2); imagesc(mean(OutAll2(:,:,p_SSD<.05),3)-mean(OutAll2(:,:,p_SSD>.05),3))
    subplot(1,3,3); imagesc(mean(AllAll2(:,:,p_SSD<.05),3)-mean(AllAll2(:,:,p_SSD>.05),3))
    suptitle(num2str(idir))
end    
% AllAll = AllIn+AllOut;
% AllAll2 = AllAll./repmat(max(max(AllAll,[],2),[],1),[size(AllAll,1) size(AllAll,2) 1]);
% InAll2 = AllIn./repmat(max(max(AllIn,[],2),[],1),[size(AllAll,1) size(AllAll,2) 1]);
% OutAll2 = AllOut./repmat(max(max(AllOut,[],2),[],1),[size(AllAll,1) size(AllAll,2) 1]);
% figure; imagesc(mean(AllAll2(:,:,smi(:,1)==1 & pSD<.05),3)-mean(AllAll2(:,:,smi(:,1)==1 & pSD>.05),3))
% figure; imagesc(mean(OutAll2(:,:,smi(:,1)==1 & pSD<.05),3)-mean(OutAll2(:,:,smi(:,1)==1 & pSD>.05),3))


