function [modu,p] = forcell_get_sig_modu_event(icell,modind,baseind,PFCreplayspikes_list)

PFCspk = PFCreplayspikes_list(PFCreplayspikes_list(:,2)==icell,:);

if size(PFCspk,1)>10

    I = PFCspk(:,3);
    nS = 500;
    ModShuff = NaN(max(I),nS);
     ModReal = NaN(max(I),1); BaseReal = ModReal;
    for i= 1:max(I)
        if sum(I==i)>0
            k = (rand(1,nS)-.5);
            j = repmat(PFCspk(I==i,1),[1 nS])+ones(sum(I==i),1)*k;
            j(j<-.5) = 1+j(j<-.5); j(j>.5) = 1-j(j>.5);            
        %         ModShuff(icell,ceil(i/2),:) = sum(j(PFCspk(I==i,2)==p(icell),:)>0 & j(PFCspk(I==i,2)==p(icell),:)<=modwin,1);
            ModShuff(i,:) = sum(j>modind(1) & j<=modind(2),1);
            q = PFCspk(I==i,1);
            ModReal(i) = sum(q>modind(1) & q<=modind(2))./range(modind);
            BaseReal(i) = sum(q>=baseind(1) & q<baseind(2))/range(baseind);
        end
    end



    meanShuf = nanmean(ModShuff);
    meanReal = nanmean(ModReal);
    Rdiff = (meanReal-meanShuf).^2;
    Sdiff = (nanmean(meanShuf)-meanShuf).^2;

    p = (sum(Sdiff>nanmean(Rdiff))+1)/(nS+1);
    modu = (nanmean(ModReal)-nanmean(BaseReal))./(nanmean(ModReal)+nanmean(BaseReal));
else p = NaN; modu = NaN;
end