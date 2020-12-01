function helper_saveandclosefig(label)

saveas(gcf,[label '.fig'],'fig')
saveas(gcf,[label '.eps'],'epsc')
saveas(gcf,[label '.tif'],'tiff')
close gcf