lowfr_pfc = 10; % Hz
highfr_pfc = 1000; % Hz
hp = 50; % Hz

wind = .040;
  shift = 20; %shift of PFC and HP spikes  
    
dt = 1/1000; % ms
nBins = 125; % 10 ms spike train
nTrials = 20; % number of simulations


hp_spikeMat1 = rand(nTrials , nBins*2 ) < hp * dt ;
hp_spikeMat = hp_spikeMat1(:,nBins+1:end);
% LowPFC_spikeMat = rand(nTrials , nBins ) < lowfr_pfc * dt ;
% LowPFC_spikeMat = [rand(nTrials , 20 ) < lowfr_pfc * dt  hp_spikeMat(:,1:105)];
LowPFC_spikeMat = hp_spikeMat1(:,[nBins+1:end]-shift);


thetacycle = false(nTrials,nBins);
thetacycle_short = thetacycle; thetacycle_short(:,1:round(nBins/2.5)) = true;
thetacycle_half = thetacycle; thetacycle_half(:,1:round(nBins/2)) = true;

figure; hold on;

subplot(2,3,4)
imagesc(LowPFC_spikeMat)
title('PFC Neuron (20ms shift from HP)')
ylabel('Simulated trial')
xlabel('Timebin (1 ms each)')
set(gca,'FontSize',18)

% subplot(3,2,2)
% imagesc(thetacycle_short)
% title('Theta cycle halves')
% ylabel('Simulated trial')
% xlabel('Timebin (1 ms each)')

subplot(2,3,2)
LowPFC_spikeMat2 = LowPFC_spikeMat; LowPFC_spikeMat2(thetacycle_short==0) = false;
imagesc(LowPFC_spikeMat2)
title('PFC, 1st half spikes')
ylabel('Simulated trial')
xlabel('Timebin (1 ms each)')
set(gca,'FontSize',18)

subplot(2,3,5)
LowPFC_spikeMat3 = LowPFC_spikeMat; LowPFC_spikeMat3(thetacycle_short==1) = false;
imagesc(LowPFC_spikeMat3)
title('PFC, 2nd half spikes')
ylabel('Simulated trial')
xlabel('Timebin (1 ms each)')
set(gca,'FontSize',18)

subplot(2,3,1)
imagesc(hp_spikeMat)
title('HP Neuron - Poisson spiking')
ylabel('Simulated trial')
xlabel('Timebin (1 ms each)')
set(gca,'FontSize',18)

a = LowPFC_spikeMat'; a = a(:);
b = hp_spikeMat'; b = b(:);
a2 = LowPFC_spikeMat2'; a2=a2(:);
a3 = LowPFC_spikeMat3'; a3=a3(:);

[cc1,lags1] = xcov(a,b,round(wind/dt),'coeff');
[cc2,lags2] = xcov(a2,b,round(wind/dt),'coeff');
[cc3,lags3] = xcov(a3,b,round(wind/dt),'coeff');
        
subplot(2,3,[3 6]); hold on
plot(lags3,cc3,'LineWidth',3)
plot(lags2,cc2,'LineWidth',3)
xlabel('Time (ms)')
ylabel('Cross-Covariance')
set(gcf,'Position', [177         268        1531         674])
set(gca,'FontSize',18)
legend('Late phase','Early phase','Location','northwest')

helper_saveandclosefig(['F:\XY_matdata\Figures\ForPaperReviews\AllCells_AllArmEvents\Theta\Simulation_PhaseDuration'])