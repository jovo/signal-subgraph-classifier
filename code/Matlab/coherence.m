fig=figure(1); clf
N=10;
blank=zeros(N);
weights=poissrnd(5,N-1,1);
fs=18;

subplot(131)
coherent=blank;
coherentedges=[1 5 9];
coherent(coherentedges(1:3),coherentedges(1:3))=reshape(weights(1:9),3,3);
imagesc(coherent)
title('coherent','FontSize',fs)
set(gca,'XTickLabel',[],'YTickLabel',[])

subplot(132)
semicoherent=blank;
semicoherentedges=[2 4 8];
semicoherent(semicoherentedges(1:3),semicoherentedges(1:3))=reshape(weights(1:9),3,3);
semicoherent(2,9)=weights(1);
semicoherent(10,8)=weights(1);
semicoherent(2,8)=0;
semicoherent(8,4)=0;
imagesc(semicoherent)
title('semicoherent','FontSize',fs)
set(gca,'XTickLabel',[],'YTickLabel',[])

subplot(133)
incoherent=blank;
edges=randperm(N^2);
incoherent(edges(1:N-1))=weights;
imagesc(incoherent)
title('incoherent','FontSize',fs)
set(gca,'FontSize',fs)

colormap(gray)


wh=[10 3];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
V.name_fig='../figs/coherence';
print('-depsc',V.name_fig)
print('-dpdf',V.name_fig)
saveas(fig,V.name_fig)