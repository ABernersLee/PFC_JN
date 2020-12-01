function [armmodu,armmodudur,tot] = make_durationmodu(thisdir,label,corrcut,jumpcut)

load(thisdir,[label '_PFCreplayspikes_binned'],[label '_replay_stnd'],[label '_replay_replayarm'],'other_cells')
load(thisdir,[label  '_replay_maxjump'],[label '_replay_corr'])

if ~exist([label '_PFCreplayspikes_binned'],'var')
    make_PFC_replayeventtriggeredmat(thisdir,label)
    load(thisdir,[label '_PFCreplayspikes_binned'])
end
eval(['PFCreplay = ' label '_PFCreplayspikes_binned;'])
eval(['maxJD = ' label '_replay_maxjump;'])
eval(['Rcorr = ' label '_replay_corr;'])
eval(['ArmM = ' label '_replay_replayarm;'])
eval(['stnd = ' label '_replay_stnd;'])
clear([label '_PFCreplayspikes_list'])

pfc = other_cells;
clear other_cells
dur = abs(diff(stnd'));

binsize = .002;
window = [-.5 .5];
ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];
modind = ind>0 & ind<=.3;
baseind = ind>=-.3 & ind<0;

armmodu = NaN(length(pfc),3);
armmodudur = armmodu;
tot = NaN(length(pfc),1);
for icell = 1:length(pfc)
    for iarm = 1:3    
        armind = ArmM==iarm;
        if sum(armind & abs(Rcorr)>corrcut & maxJD<jumpcut)>0
            m = sum(PFCreplay(icell,modind,armind & abs(Rcorr)>corrcut & maxJD<jumpcut),2);
            b = sum(PFCreplay(icell,baseind,armind & abs(Rcorr)>corrcut & maxJD<jumpcut),2); %./2;        
            rm = squeeze((m-b)./(m+b));            
            armmodudur(icell,iarm) = corr(rm,dur(armind & abs(Rcorr)>corrcut & maxJD<jumpcut)','rows','complete');
            armmodu(icell,iarm) = (nanmean(m)-nanmean(b))./(nanmean(m)+nanmean(b));
        end        
    end  
    m = sum(PFCreplay(icell,modind,abs(Rcorr)>corrcut & maxJD<jumpcut),2);
    b = sum(PFCreplay(icell,baseind,abs(Rcorr)>corrcut & maxJD<jumpcut),2); %./2;       
    rm = squeeze((m-b)./(m+b));
    tot(icell,:) = corr(rm,dur(abs(Rcorr)>corrcut & maxJD<jumpcut)','rows','complete');
end
    

    
    


