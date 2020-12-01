d2 = dir('*mat');
newcol = [75 0 130;34 139 34;255 130 0]/255;
for id = 1:size(d2,1)
   load(d2(id).name,'other_cells','InFR','OutFR','armposindex')
   for icell = 1:length(other_cells)       
       figure; hold on
       clear f g
       for iarm = 1:3
%            if iarm ==1; iarmplot = [1 2]; else; iarmplot = (iarm-1)*2+1;end
           ind = find(armposindex(:,iarm))'-find(armposindex(:,iarm),1,'first')+1;
           dat = InFR(other_cells(icell),armposindex(:,iarm));
           f(iarm) = subplot(3,2,1+(iarm-1)*2);
            patch([ind ind(end:-1:1)],[dat min(dat)*ones(size(ind))],'red','FaceColor',newcol(iarm,:),'EdgeColor',newcol(iarm,:))
            axis tight
            
            dat = OutFR(other_cells(icell),armposindex(:,iarm));
            g(iarm) = subplot(3,2,2+(iarm-1)*2);
            patch([ind ind(end:-1:1)],[dat min(dat)*ones(size(ind))],'red','FaceColor',newcol(iarm,:),'EdgeColor',newcol(iarm,:))
            axis tight
       end
       linkaxes(f,'x')
       linkaxes(g,'x')
       set(gcf,'Position',[669   672   566   217])
       set(gcf,'renderer','Painters')
       helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Basics\' d2(id).name(1:end-4) '_Cell' num2str(other_cells(icell))])
   end
end
   