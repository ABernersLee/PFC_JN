function make_OpenPlaceFields(dirs)
cd(dirs.homedir)
d2 = dir('*.mat');
binsize = 20;
SigmaField = 2;
Two_D_Filter=fspecial('gaussian',[10 10],SigmaField);
velcutoff = 5;

for id = 1:size(d2,1)
    clear OpenFR OpenFRsmoothed
    load(d2(id).name,'spikedata','vel','pos')
                        
    pos2 = pos(:,2:3);
    pos2(:,1) = ceil(pos2(:,1)/binsize);
    pos2(:,2) = ceil(pos2(:,2)/binsize);
    pos2 = pos2-min(pos2)+1;
    tm = diff(pos(:,1)); tm = [tm;tm(end)];

    pos3 = sub2ind(max(pos2),pos2(:,1),pos2(:,2));
    pos3(vel<velcutoff) = NaN;
    ind = 1:max(pos2(:,1))*max(pos2(:,2));
    occ = NaN(size(ind,2),1);
    f = zeros(size(ind,2),max(spikedata(:,2)));
    for j = 1:length(ind)
        occ(ind(j)) = sum(tm(pos3==ind(j)));
        [h,k] = histc(spikedata(pos3(spikedata(:,3)) == ind(j),2),1:max(spikedata(:,2)));
        if ~isempty(k)
            f(ind(j),k) = h(k);            
        end
        if rem(j,1000)==0
            disp(['j = ' num2str(j)])
        end
    end
    occ2 = reshape(occ,[max(pos2)]);        
    f2 = reshape(f,[max(pos2) size(f,2)]);

    OpenFR = f2./repmat(occ2,[1 1 size(f2,3)]);
    OpenFRsmoothed = NaN(size(OpenFR));
        for icell = 1:max(spikedata(:,2))            
            dat = OpenFR(:,:,icell);
            datnan = isnan(dat);
            dat(datnan) = 0;
            dat = filter2(Two_D_Filter,dat);            
            dat(datnan) = NaN;     
            OpenFRsmoothed(:,:,icell) = dat;
        end
        
    save(d2(id).name,'OpenFR','OpenFRsmoothed','-append')            
end
    