function plot_sig_subgraph(params,alg)
%% plot params

figure(1), clf
nrows=1;
ncols=3;
fs=12;


subplot(nrows,ncols,1)
imagesc(params.tru)
axis('square')
title('true','fontsize',fs)
xlabel('vertices','fontsize',fs)
ylabel('vertices','fontsize',fs)

subplot(nrows,ncols,2)
imagesc(params.coh)
axis('square')
title('coherent','fontsize',fs)


subplot(nrows,ncols,3)
imagesc(params.inc)
axis('square')
title('incoherent','fontsize',fs)



if alg.save
    wh=[5 2];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname];
    print('-dpdf',figname)
    print('-dpng',figname)
    saveas(gcf,figname)
end

