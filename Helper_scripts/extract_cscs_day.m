function extract_cscs_day(dirs,params)

% tic

cd(dirs.spikedatadir)

load([dirs.spikedatadir '\' params.ident '.mat'],'rawspikedata')    
TTs = unique(rawspikedata(:,3));
clear rawspikedata

cd(params.daydir)

st = min(min(params.Run_Times)); 
nd = max(max(params.Run_Times)); 

for itt = 1:length(TTs)

    tt = TTs(itt);
    cscname=[dir(['CSC' num2str(4*tt-3) '.ncs']);dir(['CSC' num2str(4*tt-2) '.ncs']);dir(['CSC' num2str(4*tt-1) '.ncs']);dir(['CSC' num2str(4*tt) '.ncs'])];
    SNR = NaN(size(cscname,1),1);

    for icsc = 1:size(cscname,1)

        filename = cscname(icsc).name;  

        if exist(filename,'file') && cscname(icsc).bytes~=16384 && cscname(icsc).bytes<120000000
            eval(['lfpname' num2str(icsc) '=filename;'])
            
            
            %load times, make sure they are in order, get run times
            raw_lfptimes = Nlx2MatCSC(filename,[1 0 0 0 0],0,1)/1e6;
            [raw_lfptimes,ord] = sort(raw_lfptimes);
            ind2 = raw_lfptimes>=st & raw_lfptimes<=nd;
            lfptimesA = raw_lfptimes(ind2);
            clear raw_lfptimes

            %extrapolate timepoints neuralynx didnt save
            lfpdiff = diff(lfptimesA); 
            lfpdiff = [lfpdiff lfpdiff(end)];
            lfptimesB = [0:511]'*lfpdiff/512+ones(512,1)*lfptimesA;                
            lfptimes = lfptimesB(:);                
            clear lfptimesA lfptimesB lfpdiff
            eval(['lfptimes' num2str(icsc) '= lfptimes;'])                


            %load raw samples, order, get run times
            raw_lfpsamples = Nlx2MatCSC(filename,[0 0 0 0 1],0,1);            
            raw_lfpsamplesA = raw_lfpsamples(:,ord);
            raw_lfpsamplesB = raw_lfpsamplesA(:,ind2);
            lfpsamples = raw_lfpsamplesB(:);                
            clear raw_lfpsamples raw_lfpsamplesA raw_lfpsamplesB
            
            
            eval(['lfpsamples' num2str(icsc) '= lfpsamples;'])                

            %load sampling frequencey
            Fs =Nlx2MatCSC(filename,[0 0 1 0 0],0,3,1);

            %get the theta signal to noise ratio
            L=length(lfptimes);
            NFFT=2^nextpow2(L);
            Y=fft(lfpsamples,NFFT);
            f=Fs/2*linspace(0,1,NFFT/2);
            SNR(icsc)=sum(abs(Y(f>4 & f<12,:)))./sum(abs(Y(f>59.9 & f<60.1,:)));   
            clear L T f NFFT lfpsamples lfptimes
        end
    end
        
    if sum(~isnan(SNR))>0
        [Theta_SNR,I] = max(SNR);        

        eval(['csctimes = lfptimes' num2str(I) ';'])
        eval(['cscsamples = lfpsamples' num2str(I) ';'])
        eval(['cscname = lfpname' num2str(I) ';'])
        clear lfptimes lfpsamples I 

        for icsc = 1:size(cscname,1)
            clear(['lfpname' num2str(icsc)],['lfpsamples' num2str(icsc)],['lfptimes' num2str(icsc)])
        end

        %extract the run times
        st2 = (repmat(csctimes,[1 size(params.Run_Times,1)]))'>=repmat(params.Run_Times(:,1),[1 length(csctimes)]);
        nd2 = (repmat(csctimes,[1 size(params.Run_Times,1)]))'<=repmat(params.Run_Times(:,2),[1 length(csctimes)]);
        rundat = (st2&nd2)';
        clear st2 nd2

        csctimes(sum(rundat,2)==0,:) = [];
        cscsamples(sum(rundat,2)==0,:) = [];

        save([dirs.cscdatadir '\' params.ident '_csc_tt' num2str(tt) '.mat'],'csctimes','cscsamples','cscname','Fs','Theta_SNR')    

        disp(['Done with tt ' num2str(tt)])
    end
end
disp(['Done with extract_cscs_day'])
        
% t = toc
    
% speed of indexing CSCs in different ways:
% 
% 39.6261 - adding to non indexed
% 37.4609 - making new variables with eval
% 42.8827 - indexing on first icsc