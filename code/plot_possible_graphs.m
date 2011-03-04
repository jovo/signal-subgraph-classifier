clear
fname='num_of_graphs';

%%
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


%%
figure(1), clf, hold all, fs=12;
hh(1)=subplot(131);
plot(log10(s),'color','k','linewidth',2);
set(gca,'YTick',[0:100:300],'YTickLabel',10.^[0:100:300])
ylabel('number of graphs','fontsize',fs)
xlabel('number of vertices','fontsize',fs)
title('unconstrained')

hh(2)=subplot(132); hold all
plot(log10(h),'LineStyle','-','linewidth',2,'color','k');
set(gca,'YTick',[0:100:300],'YTickLabel',10.^[0:100:300])
xlabel('number of vertices','fontsize',fs)
title('fixed size')
set(gca,'XLim',[0 100])
box on
linkaxes(hh)

subplot(133), qmax=100; 
plot(minn:qmax,log10(numcan(minn:qmax,:)),'k','linewidth',2)
xlabel('number of vertices','fontsize',fs)
set(gca,'YTick',[0:5:20],'YTickLabel',10.^[0:5:20])
title('a single star-vertex')

print_fig(['../figs/' fname],[6 2]*2)

