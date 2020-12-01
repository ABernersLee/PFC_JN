function extract_spikes_from_xclust_day(dirs,params)
cd(params.daydir)
numtetrodes = 40;
celldat = [];
CN = 0;

if strcmp(params.Rat_Name,'Chesapeake') || strcmp(params.Rat_Name,'Haw')
    disp('error, this should only be used on XY data')
end
cd('Converted')

%In Converted Folder, get tt folder names


for it = 1:numtetrodes            
    if isfolder(['tt' num2str(it)])
        cd(['tt' num2str(it)])
        d = dir;
        d2 = extractfield(d,'name');   

        %in tt folder, get cluster names
        celldr =d2(~cellfun('isempty',strfind(d2,'cl-','ForceCellOutput',1)));

        % dont load photos, mat files, or files labed as bad or
        % suspiciously long
        touse = ~contains(celldr,'bad') & ...
                ~contains(celldr,'Bad') & ...
                ~contains(celldr,'BAD') & ...
                ~contains(celldr,'.PNG') & ...
                ~contains(celldr,'Corrected') & ...
                ~contains(celldr,'.mat') & ...
                ~contains(celldr,'old') & ...
                ~contains(celldr,'missing');

        celldr = celldr(touse)
        for ic = 1:length(celldr)
            %load each cluster
            tempdat2 = load(celldr{ic});  
            if ~isempty(tempdat2)
                CN = CN+1;
                %add spike time, cell number, tetrode number, and spike widths
                tempdat = [tempdat2(:,8) CN*ones(size(tempdat2,1),1) it*ones(size(tempdat2,1),1) tempdat2(:,6)]; 
                clear tempdat2
                %add to cell data matrix
                celldat = cat(1,celldat,tempdat);
                clear tempdat
            end
        end   
        cd ../
%             disp(['Done with tt ' num2str(it)])
    end
end


%sort by time
[~,ind] = sort(celldat(:,1));
celld = celldat(ind,:);
clear celldat

%extract the run times
st = (repmat(celld(:,1),[1 size(params.Run_Times,1)]))'>=repmat(params.Run_Times(:,1),[1 length(celld(:,1))]);
nd = (repmat(celld(:,1),[1 size(params.Run_Times,1)]))'<=repmat(params.Run_Times(:,2),[1 length(celld(:,1))]);
rundat = (st&nd)';
clear st nd

rawspikedata = celld;
rawspikedata(sum(rundat,2)==0,:) = [];


save([dirs.spikedatadir '\' params.ident '.mat'],'rawspikedata','-append')
disp('Done with extract_spikes_from_xclust_day')