function make_PFC_cellwidth(thisdir)

load(thisdir, 'other_cells','rawspikedata')

PFCwidth = NaN(length(other_cells),1);
for icell = 1:length(other_cells)
    PFCwidth(icell) = mean((rawspikedata(ismember(rawspikedata(:,2),other_cells(icell)),4)./32556)*1000);
end

save(thisdir, 'PFCwidth','-append')