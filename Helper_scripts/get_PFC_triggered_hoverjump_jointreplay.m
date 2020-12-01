function get_PFC_triggered_hoverjump_jointreplay(thisdir,label)

load(thisdir,'other_cells','spikedata')
p = other_cells; clear other_cells
PFCspk = spikedata(ismember(spikedata(:,2),p),1:2); clear spikedata

load(thisdir,[label  '_seqisjumpbin'],[label  '_seqjumpdist'],[label  '_seqstepst'],[label  '_seqtimes_btwn'],[label  '_seqindex_btwn'])
eval(['isjumpbin = ' label '_seqisjumpbin;'])
eval(['stepst = ' label '_seqstepst;'])
eval(['jumpdist = ' label '_seqjumpdist;'])
eval(['seqtimes = ' label '_seqtimes_btwn;'])
eval(['seqindex = ' label '_seqindex_btwn;'])
clear([label '_seqindex_btwn'],[label '_seqisjumpbin'] , [label '_seqstepst'] ...
    , [label '_seqjumpdist'], [label '_seqtimes_btwn'])
