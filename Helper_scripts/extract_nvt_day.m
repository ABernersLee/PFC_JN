function params = extract_nvt_day(dirs,params)
cd(params.daydir)

videostring ='VT1.nvt';

%extract raw nvt file
if ismac
    eval(sprintf('[ptimes,Xpos,Ypos,HD]=Nlx2MatVT_v3(''%s'',[1 1 1 1 0 0],0,1);',videostring));    
elseif ispc
    eval(sprintf('[ptimes,Xpos,Ypos,HD]=Nlx2MatVT(''%s'',[1 1 1 1 0 0],0,1);',videostring));
end
clear videostring

%convert times to seconds
ptimes = ptimes/1e6;

%concatonate into position data matrix
pos1=[ptimes',Xpos',Ypos',HD'];
clear ptimes Xpos Ypos HD

%if times are not in order, sort them
[~,index]=sortrows(pos1(:,1));
if sum(diff(index)<0)>0        
    pos1 = pos1(index,:);
    disp('There are position times out of order')    
end
clear index



%extract the run times
st = (repmat(pos1(:,1),[1 size(params.Run_Times,1)]))'>=repmat(params.Run_Times(:,1),[1 length(pos1(:,1))]);
nd = (repmat(pos1(:,1),[1 size(params.Run_Times,1)]))'<=repmat(params.Run_Times(:,2),[1 length(pos1(:,1))]);
rundat = (st&nd)';
clear st nd

rawpos = pos1;
rawpos(sum(rundat,2)==0,:) = [];
rundat(sum(rundat,2)==0,:) = [];

params.ident = [params.Rat_Name '_' num2str(params.Date)];
save([dirs.spikedatadir '\' params.ident '.mat'],'rawpos','rundat','params','-append')

disp('Done with extract_nvt_day')
    