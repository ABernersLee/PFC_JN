eacharm1 = NaN(3,3,11); aall = NaN(3,11);
for id = 1:11
    
%     load(d2(id).name,'RP_Arm_Classifier_type1_eacharm','RP_Arm_Classifier_type1_a','RP_Arm_Classifier_type1_s')
    load(d2(id).name,'RP_Arm_ds_Classifier_type1_eacharm','RP_Arm_ds_Classifier_type1_a','RP_Arm_ds_Classifier_type1_s')
    eacharm = RP_Arm_ds_Classifier_type1_eacharm;
    a = RP_Arm_ds_Classifier_type1_a;
    figure; hold on
    subplot(1,2,1);
    imagesc(eacharm./sum(eacharm,2));hold on; axis xy; ylabel('Real'); xlabel('Guess'); c = colorbar;
    c.Label.String = 'Probability';
    subplot(1,2,2); hold on;
    bar([sum(eacharm,2) sum(eacharm)'])
    legend('Real','Guesses')
    set(gca,'xtick',[1:3])
    ylabel('Number of events')
    xlabel('Arm')
    set(gcf,'Position',[2002         369         892         291])
    helper_saveandclosefig(['F:\XY_matdata\Figures\ForPaperJN\AllCells_AllArmEvents\Classify\Sessions\ByArm_propnum_' d2(id).name(1:end-4) '_ds'])
    
    eacharm1(:,:,id) = eacharm;
    aall(:,id) = a;
end


%%
figure; histogram(squeeze(min(sum(eacharm1)./sum(sum(eacharm1)))),20)
mean(squeeze(min(sum(eacharm1)./sum(sum(eacharm1)))))


figure; hold on;
eye2 = repmat(eye(3,3)==1,[1 1 11]);
ee = eacharm1./sum(eacharm1,2);
eee = ee(eye2);
e1 = reshape(eee,[3 11]);
e1n1 = e1; e1n2 = e1;
e1n1(e1>(1/3)) = NaN;
e1n2(e1<=(1/3)) = NaN;
figure; hold on;
plot(e1n1,'ok')
plot(e1n2,'or')
ef = errorbar([1:3]+.1,mean(e1,2),std(e1,[],2)./sqrt(size(e1,2)),'k');
