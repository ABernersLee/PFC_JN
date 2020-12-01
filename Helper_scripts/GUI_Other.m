function GUI_Other(thisdir,dat)

k = figure; hold on; 
fax =gca;
set(k,'Position',[0,0,1355,995]);

% load(thisdir,'pos','spikedata')
load([thisdir '\Position_Data.mat'],'Position_Data')
load([thisdir '\Spike_Data.mat'],'Spike_Data')
pos = Position_Data; clear Position_Data;
[~,~,i1] = histcounts(Spike_Data(:,1),pos(:,1));
spikedata = [Spike_Data i1]; clear Spike_Data i1

dat.st = 1;
dat.now = 9;
dat.cellnum = 17;
dat.movebin = 1;
dat.waitbin = 0;
dat.pos = pos;
xlim([min(pos(:,2)) max(pos(:,2))])
ylim([min(pos(:,3)) max(pos(:,3))])
plot(pos(:,2),pos(:,3),'.','Color',[.5 .5 .5])
plot(pos(dat.st:dat.now-8,2),pos(dat.st:dat.now-8,3),'.k','MarkerSize',10)
plot(pos(dat.now-5:dat.now,2),pos(dat.now-5:dat.now,3),'+k','MarkerSize',20)
text(.1,-30,num2str(dat.now))


% Cell Number
u1=uicontrol;
set(u1,'Position',[10 880 100 30]);
set(u1,'Style','Edit','Tag','edit1');
set(u1,'String','93')
axes('Units','Pixels','Position',[0 900 100 30],'Visible','off');
text(.1,1,'PFC Cell:','FontSize',12)

% switch cell
u2=uicontrol;
set(u2,'Position',[10 840 100 30]);
set(u2,'UserData',dat.cellnum);
% set(u2,'Callback',@switch_cellto)
set(u2,'String','Switch Cell','FontSize',14);

%run
u3=uicontrol;
set(u3,'Position',[10 460 100 30]);
set(u3,'UserData',dat)
set(u3,'Callback',@GUIrun);
set(u3,'String','Run','FontSize',10);

%pause
u3=uicontrol;
set(u3,'Position',[10 420 100 30]);
% set(u3,'Callback',@GUIpause);
set(u3,'Value',false);
set(u3,'String','Pause','FontSize',10);
set(u3,'Style', 'togglebutton')

%waitbin
u4=uicontrol;
set(u4,'Position',[10 320 100 30]);
set(u4,'Style','Edit','Tag','editwait');
set(u4,'String','0')
text(.1,-18,'waitbin:','FontSize',12)

%movebin
u5=uicontrol;
set(u5,'Position',[10 260 100 30]);
set(u5,'Style','Edit','Tag','editmove');
set(u5,'String','1')
text(.1,-20,'movebin:','FontSize',12)

%movebin
u6=uicontrol;
set(u6,'Position',[10 380 100 30]);
set(u6,'Style','Edit','Tag','currbin');
set(u6,'String','9')



    function GUIrun(~,~)
        axes(fax)
        while true 
            dat.now = str2num(get(u6,'String'));              
            dat.pfccell = str2num(get(u1,'String'));              
            ind = spikedata(spikedata(:,2)==dat.pfccell,3);
            ind(ind>dat.now | ind<dat.st) = [];
            pause(dat.waitbin)
            cla(fax)            
            plot(pos(:,2),pos(:,3),'.','Color',[.5 .5 .5])
            plot(pos(dat.st:dat.now-5,2),pos(dat.st:dat.now-5,3),'.k','MarkerSize',10)
            plot(pos(dat.now-8:dat.now,2),pos(dat.now-8:dat.now,3),'*k','MarkerSize',20)
            plot(pos(ind,2),pos(ind,3),'or','LineWidth',2)
            dat.now = dat.now+dat.movebin;   
            dat.st = max([1;dat.now-100]);
            drawnow
            tostop = get(u3, 'Value');
            dat.movebin = str2num(get(u5,'String'));  
            dat.waitbin = str2num(get(u4,'String'));  
            
            set(u6,'String',num2str(dat.now))
            
          if tostop
            break;
          end  
        end
    end
    
    

end
