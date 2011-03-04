function plot_model(params,alg)
%% plot params


figure(1), clf
nrows=1;
ncols=3;
fs=12;

emax=max(max([params.E0(:) params.E1(:) abs(params.E0(:)-params.E1(:))]));
subplot(nrows,ncols,1), cla
image(60*(params.E0/emax))
axis('square')
colormap('gray')
title('class 0')
xlabel('vertices','fontsize',fs)
ylabel('vertices','fontsize',fs)

subplot(nrows,ncols,2), cla
image(((60*params.E1/emax)))
axis('square')
title('class 1','fontsize',fs)

subplot(nrows,ncols,3), cla
image((60*abs(params.E0-params.E1)/emax))
axis('square')
title('difference','fontsize',fs)



if alg.save
    wh=[5 2];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_params'];
    print('-dpdf',figname)
    print('-dpng',figname)
    saveas(gcf,figname)
end

