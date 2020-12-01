function helper_savefig(label)

saveas(gcf,[label '.fig'],'fig')
saveas(gcf,[label '.tif'],'tiff')
saveas(gcf,[label '.eps'],'epsc')
