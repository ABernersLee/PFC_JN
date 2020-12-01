cd('E:\XY_matdata\AllDays')
std_cutoff = 2;
d2 = dir('*.mat');
velcutoff = 5;
smoothsize = 30; withmodpatch = 1; 
for ilabel = 2 %:2
    if ilabel == 1; label = 'RP'; elseif ilabel==2; label = 'SD'; end
    
%     if 0
    for id =  4:11
%         try
        thisdir = d2(id).name;
% 
        if id>4
        make_directionalfields(thisdir,velcutoff)
        disp('1')
        
%             
        if ilabel == 2
            make_SDCandEvents_spikedensity(thisdir,std_cutoff) 
        elseif ilabel == 1
            make_RPCandEvents_ripple(thisdir,std_cutoff) 
        end
        disp('1')    
      
        make_PosteriorfromCandEvents(thisdir,label)        
      
        disp('2')
        
      
        make_singleandjointreplay_fromPosterior(thisdir,label,5000) %added shuffle
        disp('3')    


    %         
        make_PFC_candeventtriggeredmat(thisdir,label);
        disp('4')            
        make_PFC_replayeventtriggeredmat(thisdir,label) % changed to be significant events only
        disp('5')    
    %         
        
        get_PFC_armtriggered_modusig(thisdir,label); %changed to be significant events only
        disp('6')    
        end
% 
%         get_PFC_fwdrevtriggered_modusig(thisdir,label)
%         disp('8')    
%         plot_PFC_Event_triggered(thisdir,label,withmodpatch,smoothsize)
% %         disp('7')
% %         end

        plot_PFC_ArmReplay_triggered(thisdir,label,withmodpatch,smoothsize)
%         disp('7')
% %         end
        
        run_all_linear_classifiers(thisdir,label,1) %updated with correct significance tests
        disp('8')
        
        plot_PFC_JointReplay_triggered(thisdir,label,withmodpatch,smoothsize,timetrig)

%         load(thisdir,'spikedata','pos','HP_Theta')        
%         HPTheta_spikephase = make_HPtheta_spikephase2(spikedata,HP_Theta,pos);
%         save(thisdir,'HPTheta_spikephase','-append')
%         disp('13')
% 
%         make_PFC_behavechange_eventtriggeredmat(thisdir,1) 
%         disp('6/10')
%         
%         prospective_coding(thisdir)
%         disp('7/10')        
        if id~=2
            for ilaptype = 3 %1:3
               get_lapbylap_PFC(thisdir,label,ilaptype,10)
%                disp([num2str(ilaptype)])
            end                
        end
%         
        disp(['Done with Day ' num2str(id)])      
%         catch ME
%             disp(['ID: ' ME.identifier])    
%             msgString = getReport(ME);
%             disp(msgString)
%             disp('Error Occured')
%         end

    end
    

    %todo  
    compare_arm_replay_theta(label) % very rough, come back to


    %done - check label is in save fig
%     PFCcell_HPtheta_byarm(velcutoff,5)
%     end

%     for ilaptype = 3 %1:3
       plot_lapbylap_PFC_all(label,ilaptype)       
%     end
    proportion_replay_modulated_all(label)  
    for itype = 1
        linear_classifier_of_replay_all(label,'Arm',itype)
        linear_classifier_of_replay_all(label,'Dir',itype)
        linear_classifier_of_replay_all(label,'RatPos',itype)
    end

%     plot3PFC_arm_modulation(label) %FIX    
% %     choicepoint_allPFC(label) % need to add to
%     plot_Velocity_Vs_Value(label)
% % run_prospective_coding_fields.m

    
    num_shuffles = 1000;
    for area = 1:2
        hs = NaN(3,3,11); ps = NaN(11,3); ps2 = NaN(11,3); 
        for id = 1:11    
                [ps(id,1),ps2(id,:),ps(id,2),ps(id,3),hs(:,:,id)] = make_jointreplayeventtriggeredmat(d2(id).name,label,num_shuffles,area);      
        end
        save(['E:\XY_matdata\AllDaysMat\predictjointreplays_20181012_' label '_area' num2str(area)],'hs','ps','ps2','num_shuffles')
    end
    
end


% singlecell_jointreplay_stemdiff(thisdir,label,num_shuffles,area)
% predictjointreplays.m