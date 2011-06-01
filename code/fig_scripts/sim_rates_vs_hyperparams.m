%% synthetic
clear, clc


% # neurons
V=70;

% signal subgraph
m=1;
s=20;

Ess = randperm(70);
Ess = Ess(1:s);
thesem = zeros(V);
vstar=round(V/2);
thesem(vstar,Ess) = 1;
thesem(Ess,vstar) = 1;
thesem=tril(thesem,-1);
Ess=find(thesem);

% parameters
p=0.1;
q=0.3;
E0=p*ones(V);
E1=p*ones(V);
E1(Ess)=q;

ntrn=200;
ntst=500;

fname=['coherent_image_V' num2str(V) '_s' num2str(s) '_p' num2str(round(p*100)) '_q' num2str(round(q*100)) '_nTr' num2str(ntrn) '_nTe' num2str(ntst)];


%% generate data, classify, and get stats
clear alg out Outtemp ns h
% alg stuff
alg(1).name='coherent';
alg(1).star_list=1:V/2;
for m=1:length(alg(1).star_list)
    alg(1).edge_list{m}=1:V*alg(1).star_list(m)/2;
end

alg(2).name='incoherent';
alg(2).edge_list=1:choose(V,2);

alg(3).name='naive bayes';
alg(3).edge_list=find(ones(V)-diag(ones(V,1)));

nAlgs=numel(alg);

%% generate data

% training data
Atrn = nan(V,V,ntrn);
Atrn(:,:,1:2:ntrn)=repmat(E0,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
Atrn(:,:,2:2:ntrn)=repmat(E1,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
for nn=1:ntrn, Atrn(:,:,nn)=tril(Atrn(:,:,nn),-1); end
ytrn=repmat([0 1],1,ntrn/2);

% testing data
Atst = nan(V,V,ntst);
Atst(:,:,1:2:ntst)=repmat(E0,[1 1 ntst/2]) > rand(V,V,ntst/2);
Atst(:,:,2:2:ntst)=repmat(E1,[1 1 ntst/2]) > rand(V,V,ntst/2);
for nn=1:ntst, Atst(:,:,nn)=tril(Atst(:,:,nn),-1); end
ytst=repmat([0 1],1,ntst/2);


%% classify
Out = plugin_classifier_cv_loop(Atrn,ytrn,alg,'HoldOut',Atst,ytst);

for i=1:nAlgs
    Lhats{i}=Out(i).Lhat;
end

save(['../../data/' fname])

%% plot rates

% load(['../../data/' fname])
clear h
figure(1), clf

gray1=0.5*[1 1 1];
gray2=0.75*[1 1 1];
sa=0.5;
fs=12;
ms=24;
xmax=1000;
h(1)=subplot(211);  hold all
semilogx(Lhats{2},'color','k','linewidth',2)
min2=min(Lhats{2});
shat_inc=find(Lhats{2}==min2);
set(gca,'XScale','log'),
xlabel('log size of signal subgraph','fontsize',fs)
ylabel('misclassification rate','fontsize',fs)
title('incoherent estimator','fontsize',fs)
axis([1 xmax 0 0.5])
set(gca,'XTick',[1 shat_inc 100 1000])

% plot nb
plot(xmax,Lhats{3},'.','color',0.5*[1 1 1],'markersize',ms)
text(xmax/3,Lhats{3}*1.1,['$\hat{L}_{nb}=$' num2str(min2)],'interp','latex')

% plot Lhat_inc
plot(shat_inc,min2,'.','color',0.5*[1 1 1],'markersize',ms)
text(shat_inc*1.05,min2*0.9,['$\hat{L}_{inc}=$' num2str(min2)],'interp','latex')
plot([shat_inc shat_inc],[0 min2],'--k')
plot([1 shat_inc],[min2 min2],'--k')

% plot Lstar
Lstar=0.1;
% semilogx([1 length(Lhats{2})],[Lstar Lstar],'--k')
plot(s,Lstar,'.','color',0.5*[1 1 1],'markersize',ms)
text(s+4,Lstar,['$\hat{L}_{*}=$' num2str(Lstar)],'interp','latex')

% plot L_chance
plot(1,0.5,'.','color',0.5*[1 1 1],'markersize',ms)
text(1.1,0.48,['$\hat{L}_{chance} \quad$  $=0.5$'],'interp','latex')


h(2)=subplot(212);
cb=imagesc(Lhats{1});
min1=min(Lhats{1}(:));
coh_inds=find(Lhats{1}==min1);
[II JJ]=ind2sub(size(Lhats{1}),coh_inds);
[I Iind]=min(II);
miny=min(Lhats{1}(:));
maxy=max(Lhats{1}(:));
colorbar('YTick',[.1:.1:.5])
xlabel('size of signal subgraph','fontsize',fs)
ylabel('number of star-vertices','fontsize',fs)
tit=['coherent estimator'];
title(tit,'interp','none','fontsize',fs)
colormap('gray')

annotation('textarrow',[0.67 0.58],0.61*[1 1],'color','w')

text(600,14,['$\hat{L}_{coh}=$' num2str(min1)],'interp','latex','color','w')

set(h,'fontsize',fs)

print_fig(['../../figs/' fname],[4 6])




