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


fname=['coherent_image_V' num2str(V) '_s' num2str(s) '_p' num2str(round(p*100)) '_q' num2str(round(q*100))];


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
nTrials=20;
ntrn=100;

%% generate data

disp(['trial ', num2str(j)])
% training data
Atrn = nan(V,V,ntrn);
Atrn(:,:,1:2:ntrn)=repmat(E0,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
Atrn(:,:,2:2:ntrn)=repmat(E1,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
for nn=1:ntrn, Atrn(:,:,nn)=tril(Atrn(:,:,nn),-1); end
ytrn=repmat([0 1],1,ntrn/2);

% testing data
ntst=1000;
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

save(['../data/' fname])

%% plot rates

% load(['../data/' fname])
clear h
figure(1), clf

gray1=0.5*[1 1 1];
gray2=0.75*[1 1 1];
sa=0.5;
fs=10;
xmax=1000;
h(1)=subplot(121);  
semilogx(Lhats{2},'color','k','linewidth',2)
xlabel('log size of signal subgraph','fontsize',fs)
ylabel('misclassification rate','fontsize',fs)
title('incoherent estimator')
axis([1 1000 0 0.5])

h(2)=subplot(122);
imagesc(Lhats{1}), colorbar
xlabel('size of signal subgraph','fontsize',fs)
ylabel('number of star-vertices','fontsize',fs)
title('coherent estimator')

set(h,'fontsize',fs)
%%
print_fig(['../figs/' fname '_Lhats'],[8 3]*1.5)


