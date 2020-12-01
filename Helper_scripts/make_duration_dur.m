function rss = make_duration_dur(thisdir,label,corrcut,jumpcut)

load(thisdir,[label '_PFCreplayspikes_list'],[label '_replay_stnd'],[label '_replay_replayarm'],'other_cells')
load(thisdir,[label  '_replay_maxjump'],[label '_replay_corr'])

if ~exist([label '_PFCreplayspikes_list'],'var')
    make_PFC_replayeventtriggeredmat(thisdir,label)
    load(thisdir,[label '_PFCreplayspikes_list'])
end
eval(['PFCreplay = ' label '_PFCreplayspikes_list;'])
eval(['maxJD = ' label '_replay_maxjump;'])
eval(['Rcorr = ' label '_replay_corr;'])
eval(['ArmM = ' label '_replay_replayarm;'])
eval(['stnd = ' label '_replay_stnd;'])
clear([label '_PFCreplayspikes_list'])

PFCreplay(maxJD>jumpcut | abs(Rcorr)<corrcut,:) = [];


pfc = other_cells;
clear other_cells
dur = abs(diff(stnd'));

rss = NaN(length(pfc),1);
if ~isempty(PFCreplay)
    for icell = 1:length(pfc)
      alldat = NaN(max(PFCreplay(:,3)),2);
      for ireplay = 1:max(PFCreplay(:,3))
          spks = sum(PFCreplay(PFCreplay(:,2)==pfc(icell) & PFCreplay(:,3)==ireplay & PFCreplay(:,1)<dur(ireplay) & PFCreplay(:,1)>0));
          if spks>0
             alldat(ireplay,1) = spks;
             alldat(ireplay,2) = dur(ireplay);
          end
      end
      rss(icell,1) = corr(alldat(:,1)./alldat(:,2),alldat(:,2),'rows','complete');
      rss(icell,2) = corr(alldat(:,1),alldat(:,2),'rows','complete');
    end
end

    
    


