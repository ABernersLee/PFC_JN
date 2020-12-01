function [prec_hp,prec_pfc,sig_pfc,toplot,toplotind] = get_ArmReplay_triggered_spatialprecision(thisdir,igroup)

%this used to be longer, now make make_ArmReplay_triggered_spatialprecision
%that saves out spatial precision

load(thisdir,'other_cells','other_cells_touse','RP_pSSDarm','spatialprecision','hpinterneurons','hp_cells')
p_SSD = RP_pSSDarm;
pfc = other_cells(other_cells_touse(:,igroup));
sig_pfc = p_SSD(other_cells_touse(:,igroup));
hp = hp_cells(~ismember(hp_cells,hpinterneurons));
prec_hp = spatialprecision.prec(hp,:);
prec_pfc = spatialprecision.prec(pfc,:);
toplot = spatialprecision.toplot;
toplotind1 = spatialprecision.toplotind;

toplotind = NaN(size(toplotind1,1),4);
toplotind(hp,1) = 1;
toplotind(pfc,1) = 2;
toplotind(other_cells(other_cells_touse(:,igroup)),2) = p_SSD(other_cells_touse(:,igroup));
toplotind(:,3:4) = toplotind1;
