function print_fig(figname,wh)

set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
print('-dpdf',figname)
print('-dpng',figname)
saveas(gcf,figname)

end


