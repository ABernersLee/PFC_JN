% Goal: re-organize all the PFC data into a better format for the paper.

flow:

%run at start and never again
transferdata.m (should only run once - should never overwrite or append)
%each loops through all the days
transferLFPandTT.m (appends)

%each one loops through all the days
make_localgamma_spikephase(type) %1 for slow 2 for fast
make_HPtheta_spikephase

~~~~~~
% make overarching script that loops through all days for these
 %ready for speciefic questions for figures after this

make_directionalfields(thisdir)
make_SDCandEvents_spikedensity(thisdir) % same as Ting
make_RPCandEvents_ripple(thisdir,rip_thresh) % same as brad
make_PosteriorfromCandEvents(thisdir,label)
make_singleandjointreplay_fromPosterior(thisdir,label)
make_hoverjump(thisdir,label)

make_slowgama_ofsteps(thisdir,label)
make_PFC_cellwidth(thisdir)
make_PFC_candeventtriggeredmat(thisdir,label)
make_PFC_replayeventtriggeredmat(thisdir,label)

get_PFC_armtriggered_modusig(thisdir,label)
get_PFC_fwdrevtriggered_modusig(thisdir,label)

linear_classifer_of_replay_withPFC(thisdir,label,type)
	%type 1 = Linear Discriminant Analysis
	%type 2 = Multinomial Logistic Regression

get_PFC_hoverjump_spikeprob_mrv(thisdir,label) *not done

make_PFC_behavechange_eventtriggeredmat(thisdir)
**HERE finish and add across laps, then runsome.m

~~~~~
% plot

get_PFC_triggered_hoverjump_jointreplay(thisdir,label) *not started
get_PFC_triggered_hoverjump_(thisdir,label) * not started

*add signififance to it
plot_PFC_Event_triggered(thisdir,label,withmodpatch,smoothsize)
plot_PFC_ArmReplay_triggered(thisdir,label,withmodpatch,smoothsize)
plot_PFC_FwdRevReplay_triggered(thisdir,label,withmodpatch,smoothsize) 

plot3PFC_arm_modulation * not saved
~~~
%analysis over all days

PFCcell_HPtheta_byarm(velcutoff,binsize) * draft, has figs and spits out arm data
compare_arm_replay_theta(label) *draft

proportion_replay_modulated_all(label)
linear_classifier_of_replay_all(label,1)
linear_classifier_of_replay_all(label,2)
plot3PFC_arm_modulation(label)

~~~~~
functions:
transferdata: Takes the data that I already processed for the replay part of the paper and transfers it over here.
		This data is only up to MakePrecessedData_ABL20180304, not the replay. I will analyze the replay differently
		for this part of the paper.
transferLFPandTT: Compares cell and tetrode info between the two places I'm getting data to make sure nothing strange is
		happening, and then saves the mean Ripple, Theta and Slow Gamma HP that I had calculated before, and the TT numbers.

make_slowgamma_spikephase: Gets the local CSC for each cell and makes a matrix the same size as spikedata
		which has slow gamma power, Zscored power, and phase of each tetrode

