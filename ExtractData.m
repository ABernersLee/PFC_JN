function ExtractData(dirs)

%Place dirs outside of homedir folder and load to run ExtractData and ProcessData
% homedir = 'E:\XY_matdata\AllDays\';
% dirs.spikedatadir = homedir;
% dirs.cscdatadir = [homedir 'cscdata\'];
% dirs.paramdir = [homedir 'params\'];
% dirs.homedir = homedir;
% cd('E:\XY_matdata\')
% save([homedir 'dirs.mat'],'dirs')
% clearvars -except dirs


cd(dirs.homedir)
d2 = dir('*.mat');
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    
    load([dirs.homedir thisdir],'params')
            
    %extracts raw position and rundata, not processed yet    
    params = extract_nvt_day(dirs,params);

    if 0
        %extracts xclust spike data and saves into each run
        extract_spikes_from_xclust_day(dirs,params)

        %extract the cscs (only one per tetrode - theta/noise ratio)
        extract_cscs_day(dirs,params)    
    end
    
    disp(['Done with extraction for day ' num2str(id)])
end