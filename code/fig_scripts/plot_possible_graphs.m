clear
fname='num_of_graphs';

%% generate number of graphs under different constraints

maxn=100;
h=nan(maxn,4950);

% generate unconstrained and size constrained
for n=2:maxn;
    d(n)=nchoosek(n,2);
    s(n)=2^d(n);
    for k=min(n,10):10:min(100,d(n)/2); %round(logspace(0,log10(d(n)/2),10));
        h(n,k)=nchoosek(d(n),k);
    end
end


% generate coherent constrained
minn=3;
maxn=100;
numcan=zeros(maxn,5);
for n=minn:maxn
    for i=5:n
        for j=i+1:n
            for k=1:5
                numcan(n,k)=numcan(n)+nchoosek((i-1)*(j-1),k);
            end
        end
    end
end

save(['../data/' fname])
%% generate examples


n=10;
m=6;
p=0.2;
q=0.5;

diags=1:n+1:n^2;

M_inc=randperm(n^2);
M_inc=M_inc(1:m);

M_cross=randperm(n);
M_cross=M_cross(1:m);

ER=ones(n)-eye(n);

INC=zeros(n);
INC(M_inc)=1;
INC=INC+INC';
INC(INC>1)=1;

INC=zeros(n);
INC(M_inc)=1;
INC=INC+INC';
imagesc(INC)


CROSS=zeros(n);
CROSS(5,M_cross(1:m/2))=1;
CROSS=CROSS+CROSS';
CROSS(CROSS>1)=1;



%% generate fig

clear hh
figure(1), clf
nrows=2;
ncols=3;
fs=12;



% figure(1), clf, hold all, fs=12;
hh(1)=subplot(231);
plot(log10(s),'color','k','linewidth',2);
set(gca,'YTick',[0:100:300],'YTickLabel',10.^[0:100:300])
ylabel('number of graphs','fontsize',fs)
xlabel('number of vertices','fontsize',fs)
% title('Complete','fontsize',fs)
title('unconstrained')

hh(2)=subplot(232); hold all
plot(log10(h),'LineStyle','-','linewidth',2,'color','k');
set(gca,'YTick',[0:100:300],'YTickLabel',10.^[0:100:300])
xlabel('number of vertices','fontsize',fs)
title('incoherent','fontsize',fs)
% title('fixed size')
set(gca,'XLim',[0 100])
box on
linkaxes(hh)

hh(3)=subplot(233); qmax=100; 
plot(minn:qmax,log10(numcan(minn:qmax,:)),'k','linewidth',2)
xlabel('number of vertices','fontsize',fs)
set(gca,'YTick',[0:5:20],'YTickLabel',10.^[0:5:20])
title('coherent','fontsize',fs)
% title('a single star-vertex')


emax=q;
hh(4)=subplot(234); cla
imagesc(ER)
colormap('gray')
xlabel('vertices','fontsize',fs)
ylabel('vertices','fontsize',fs)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])
axis('square')

hh(5)=subplot(235); cla
imagesc(INC)
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])
axis('square')

% text(1.0,-1.5,'Signal Subgraph Types','fontsize',fs+4)

hh(6)=subplot(236); cla
imagesc(CROSS)
colormap('gray')
set(gca,'fontsize',fs,'DataAspectRatio',[1 1 1])
axis('square')

set(h,'FontSize',fs)

%%
print_fig(['../figs/' fname],[5 2]*2)


% 
% wh=[7.5 3.5]*2;   %width and height
% set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
% figname='~/Research/figs/misc/com_inc_coh';
% print('-dpdf',figname)
% print('-dpng',figname)
% saveas(gcf,figname)

