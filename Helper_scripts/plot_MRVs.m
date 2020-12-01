function plot_MRVs(dirs)

cd(dirs.homedir)
d2 = dir('*.mat');
mrv_mu_candsig_modu = [];
mrv_mu_hp = [];
if ~isfolder(['E:\XY_matdata\Figures\ForPaper\Theta\MRVs\'])
    mkdir(['E:\XY_matdata\Figures\ForPaper\Theta\MRVs\'])
end
plotexamples = true;
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'HPTheta_spikephase','vel','theta_globalzero','hp_cells','other_cells','hpinterneurons','RP_Cand_sig_modu_include','RP_pSSDarm')
    theta = HPTheta_spikephase;
    theta(:,6) = mod(theta(:,6)-theta_globalzero,360);
    theta(vel(theta(:,3))<5,:) = [];
    
    hp = hp_cells(~ismember(hp_cells,hpinterneurons));
    pfc = other_cells;


    mrv_mu_pfc = NaN(length(pfc),2);
    for icell = 1:length(pfc)    
        dat = deg2rad(theta(theta(:,2)==pfc(icell),6));
        mrv_mu_pfc(icell,1) = circ_r(dat);
        mrv_mu_pfc(icell,2) = abl_circ_mean_ci(dat,[],[],[],[],1);        
        
        if plotexamples
            [t,r] = rose(dat,20);
            polarplot(t,r)
            cp = circ_plot(dat,'hist',[],20,true,true,'linewidth',3,'color','r');            
            cp.Children(2).Color = 'k';
            cp.Children(2).LineWidth = 2;
            title([thisdir ' Cell ' num2str(pfc(icell))])
            helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Theta\MRVs\Example_' thisdir '_Cell' num2str(pfc(icell))])
        end
    end

    mrv_mu_candsig_modu = [mrv_mu_candsig_modu; mrv_mu_pfc RP_Cand_sig_modu_include(:,1:2)];
    
    mrv_mu_hp1 = NaN(length(hp),2);
    for icell = 1:length(hp)    
        dat = deg2rad(theta(theta(:,2)==hp(icell),6));
        mrv_mu_hp1(icell,1) = circ_r(dat);
        mrv_mu_hp1(icell,2)  =abl_circ_mean_ci(dat,[],[],[],[],1);
    end
    mrv_mu_hp = [mrv_mu_hp; mrv_mu_hp1];
end
%%
figure; hold on
dat1 = mrv_mu_candsig_modu(:,1); dat2 = mrv_mu_candsig_modu(:,4);
% plot(dat1,abs(dat2),'.')
plot(log(dat1),abs(dat2),'k.','MarkerSize',20)
[r1,p1] = corr(log(dat1),abs(dat2));
[r2,p2] = corr(dat1,abs(dat2));
b = polyfit(log(dat1),abs(dat2),1);
x1 = [min(log(dat1)) max(log(dat1))];
y2 = polyval(b,x1);
plot(x1,y2,'r','LineWidth',3)
text(-5,.6,['Log-Lin corr, r = ' num2str(round(r1,2)) ' p = ' num2str(round(p1,2,'significant'))],'Color','r')
text(-5,.55,['Linear corr, r = ' num2str(round(r2,2)) ' p = ' num2str(round(p2,2))],'Color','r')
xlabel('Log(MRV to HP Theta)')
ylabel('Ripple Modulation')
set(gca,'FontSize',18)
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Theta\MRVs\RippleVsMRV'])
%%

dat = deg2rad(mrv_mu_candsig_modu(:,2));
[t,r] = rose(dat,20);
polarplot(t,r);
cp = circ_plot(dat,'hist',[],20,true,true,'linewidth',3,'color','r');
cp.Children(2).Color = 'k';
cp.Children(2).LineWidth = 2;
set(gca,'FontSize',18)
title('All mPFC Cells')
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Theta\MRVs\AllPFC'])

dat = deg2rad(mrv_mu_candsig_modu(mrv_mu_candsig_modu(:,3)==1,2));
[t,r] = rose(dat,20);
polarplot(t,r);
cp = circ_plot(dat,'hist',[],20,true,true,'linewidth',3,'color','r');
cp.Children(2).Color = 'k';
cp.Children(2).LineWidth = 2;
set(gca,'FontSize',18)
title('All Ripple Modulated mPFC Cells')
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Theta\MRVs\AllPFC_RipMod'])

dat = deg2rad(mrv_mu_hp(:,2));
[t,r] = rose(dat,20);
polarplot(t,r);
cp = circ_plot(dat,'hist',[],20,true,true,'linewidth',3,'color','r');
cp.Children(2).Color = 'k';
cp.Children(2).LineWidth = 2;
set(gca,'FontSize',18)
title('All HP Cells')
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Theta\MRVs\AllHP'])