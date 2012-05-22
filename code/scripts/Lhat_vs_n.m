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

savestuff=0;

fname=['homo_V' num2str(V) '_s' num2str(s) '_p' num2str(round(p*100)) '_q' num2str(round(q*100))];


%% plot params and samples
figure(2), clf
fs=10;

h(2)=subplot(222);
image(E1*64+0.5)
axis('square');
op2b=get(gca,'Position');
colorbar('YTick',[0.01 min(p,q) max(p,q) 1]*64.5,'YTickLabel',[0 min(p,q) max(p,q) 1])
set(gca,'Position',op2b)
title('Class 1','fontsize',fs)

h(1)=subplot(221);
image(E0*64.5)
colormap('gray')
axis('square');
title('Class 0','fontsize',fs)
xlabel('vertex','fontsize',fs)
ylabel('expected value')
op1=get(gca,'Position');
op1(3:4)=op2b(3:4);
set(gca,'Position',op1)


h(3)=subplot(223);
imagesc(tril(E0 > rand(V),-1))
axis('square');
ylabel('sample graph')
xlabel('vertex','fontsize',fs)

h(4)=subplot(224);
imagesc(tril(E1 > rand(V),-1))
axis('square');


for i=1:numel(h)
    set(h(i),'fontsize',fs)
end

if savestuff==1
    print_fig(['../figs/' fname '_mean_sample'],[6 4])
end

%% simulate data for coherogram plots
clear h SigMat coherent incoherent coherogram Coherogram wset

ns=[4 16 64 256]; %[20 50 100 300];
ncols=length(ns);


for i=1:ncols
    
    % training data
    ntrn=ns(i);
    Atrn = nan(V,V,ntrn);
    Atrn(:,:,1:2:ntrn)=repmat(E0,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
    Atrn(:,:,2:2:ntrn)=repmat(E1,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
    for nn=1:ntrn, Atrn(:,:,nn)=tril(Atrn(:,:,nn),-1); end
    
    ytrn=repmat([0 1],1,ntrn/2);
    
    constants   = get_constants(Atrn,ytrn);     % get constants to ease classification code
    SigMat{i}   = run_get_fisher_pvals(Atrn,constants);
    incoherent{i} = incoherent_estimator(SigMat{i},s);
    [coherent{i} cocount{i}]= coherent_estimator(SigMat{i},m,s);
    [coherogram{i} Coherogram{i} wset{i}] = get_coherograms(Atrn,ytrn);
    
end


%% plot SigMat and coheregrams
% load(['../../data/results/' fname])

figure(4), clf
wh = [.21 .18];
side = 0.07;
bottom=0.07;
rowpos=0.22;
colpos=0.21;
fs=6;
nrows=4;


for i=1:ncols
    
    
    h(1+(i-1)*nrows) = subplot('Position',[0*colpos+side (4-i)*rowpos+bottom wh]);
    
    imagesc(-log10(SigMat{i}))
    axis('square')
    ylabel(['n=', num2str(ns(i))],'fontsize',fs)
    if i<ncols,
        set(gca,'Xtick',[],'YTick',[]),
    else
        xlabel('vertex','fontsize',fs)
%         ylabel([{['n=', num2str(ns(i))]}; {'vertex'}],'fontsize',fs)
        ylabel([{['n=', num2str(ns(i))]}],'fontsize',fs)
        set(gca,'YTick',[]),
    end
    if i==1, title([{'negative log'};  {'significance matrix'}],'fontsize',fs), end
    
    h(2+(i-1)*nrows) = subplot('Position',[1*colpos+side (4-i)*rowpos+bottom wh]);
    ass=zeros(V); ass(incoherent{i})=1;
    imagesc(ass)
    axis('square')
    if i<ncols, set(gca,'Xtick',[],'YTick',[]),
    else set(gca,'Ytick',[]), xlabel('vertex','fontsize',fs)
    end
    if i==1, title([{'incoherent'}; {'estimate'}],'fontsize',fs), end
    ylabel(['# correct = ', num2str(length(intersect(incoherent{i},Ess)))],'fontsize',fs)
    
    h(3+(i-1)*nrows) = subplot('Position',[2*colpos+side (4-i)*rowpos+bottom wh]);
    ass=zeros(V); ass(coherent{i})=1;
    imagesc(ass)
    if isempty(coherent{i}), image(ass); end
    axis('square')
    if i<ncols, set(gca,'Xtick',[],'YTick',[]),
    else set(gca,'Ytick',[]), xlabel('vertex','fontsize',fs)
    end
    if i==1, title([{'coherent'}; {'estimate'}],'fontsize',fs), end
    ylabel(['# correct = ', num2str(length(intersect(coherent{i},Ess)))],'fontsize',fs)
    
    h(4+(i-1)*nrows) = subplot('Position',[3*colpos+side (4-i)*rowpos+bottom wh]);
    imagesc(Coherogram{i}(:,1:cocount{i}))
    axis('square')
    set(gca,'Ytick',[])
    if i==1, title([{'coherogram'}],'fontsize',fs), end
    if i==ncols, xlabel('log significance','fontsize',fs), end
    if cocount{i}>2,
        set(gca,'XTick',ceil(linspace(2,cocount{i},3)),'XTickLabel',...
            0.1*round(10*log10((wset{i}(ceil(linspace(2,cocount{i},3))))))),
    end
    
end
colormap('gray')
for i=1:length(h), set(h(i),'fontsize',fs), end

if savestuff==1
    print_fig(['../../figs/' fname '_coherogram'],[5 4])
end



%% get Lstar

ntst=1000;
Atst = nan(V,V,ntst);
Atst(:,:,1:2:ntst)=repmat(E0,[1 1 ntst/2]) > rand(V,V,ntst/2);
Atst(:,:,2:2:ntst)=repmat(E1,[1 1 ntst/2]) > rand(V,V,ntst/2);
ytrn=repmat([0 1],1,ntst/2);

params.lnE0=log(E0);
params.lnE1=log(E1);
params.ln1E0=log(1-E0);
params.ln1E1=log(1-E1);
params.lnprior0=log(1/2);
params.lnprior1=log(1/2);

subspace=find(ones(V)-diag(ones(V,1)));

bayes_incorrects=nan(ntst,1);
for i=1:ntst
    bayes_incorrects(i) = plugin_bern_classify(Atst(:,:,i),params,subspace,ytrn(i));
end
Lstar=mean(bayes_incorrects)

%% generate data, classify, and get stats
clear alg out Outtemp ns h
% alg stuff
alg(1).name='coherent';
alg(1).star_list=m;
alg(1).edge_list{1}=s;

alg(2).name='incoherent';
alg(2).edge_list=s;

alg(3).name='naive bayes';
alg(3).edge_list=find(ones(V)-diag(ones(V,1)));

nAlgs=numel(alg);
nTrials=20;
ns=[10:10:50, 100:50:300];


for i=1:length(ns)
    
    disp(['ntrn ', num2str(ns(i))])
    
    for j=1:nTrials
        
        disp(['trial ', num2str(j)])
        % training data
        ntrn=ns(i);
        Atrn = nan(V,V,ntrn);
        Atrn(:,:,1:2:ntrn)=repmat(E0,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
        Atrn(:,:,2:2:ntrn)=repmat(E1,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
        for nn=1:ntrn, Atrn(:,:,nn)=tril(Atrn(:,:,nn),-1); end
        
        ytrn=repmat([0 1],1,ntrn/2);
        
        % testing data
        ntst=100;
        Atst = nan(V,V,ntst);
        Atst(:,:,1:2:ntst)=repmat(E0,[1 1 ntst/2]) > rand(V,V,ntst/2);
        Atst(:,:,2:2:ntst)=repmat(E1,[1 1 ntst/2]) > rand(V,V,ntst/2);
        for nn=1:ntst, Atst(:,:,nn)=tril(Atst(:,:,nn),-1); end
        ytst=repmat([0 1],1,ntst/2);
        
        
        % classify
        Outtemp = plugin_classifier_cv_loop(Atrn,ytrn,alg,'HoldOut',Atst,ytst);
        
        % get stats
        Out(1).Lhat(i,j)   = Outtemp(1).Lhat;
        Out(1).correct_edges(i,j) = length(intersect(Outtemp(1).subspace{1},Ess));
        
        Out(2).Lhat(i,j)   = Outtemp(2).Lhat;
        Out(2).correct_edges(i,j) = length(intersect(Outtemp(2).subspace{1},Ess));
        
        Out(3).Lhat(i,j)   = Outtemp(3).Lhat;
    end
    
end


for n=1:3
    keeps(n).correct_edges=Out(n).correct_edges;
    keeps(n).Lhat=Out(n).Lhat;
    
    keeps(n).misedge_rate=1-mean(keeps(n).correct_edges,2)/s;
    keeps(n).misedge_ste=std(keeps(n).correct_edges,[],2)/(s*sqrt(nTrials));
    keeps(n).Lhat_avg=mean(keeps(n).Lhat,2);
    keeps(n).Lhat_ste=std(keeps(n).Lhat,[],2)/sqrt(nTrials);
end

if savestuff==1
    save(['../../data/' fname])
end

%% plot rates
% load(['../../data/' fname])

figure(1), clf

gray1=0.5*[1 1 1];
gray2=0.75*[1 1 1];
sa=0.5;
fs=12;


%
h(1)=subplot(211); hold all
errorbar(ns,keeps(1).misedge_rate(1:length(ns)),keeps(1).misedge_ste(1:length(ns)),'k','linewidth',2)
errorbar(ns,keeps(2).misedge_rate(1:length(ns)),keeps(2).misedge_ste(1:length(ns)),'color',gray1,'linewidth',2)
axis([0 ns(end) 0 1])
% set(gca,'YTick',[2:2:s])
xlabel('# training samples','fontsize',fs)
ylabel('missed-edge rate','fontsize',fs)

h(2)=subplot(212); hold all
errorbar(ns,keeps(1).Lhat_avg(1:length(ns)),keeps(1).Lhat_ste(1:length(ns)),'k','linewidth',2)
errorbar(ns,keeps(2).Lhat_avg(1:length(ns)),keeps(2).Lhat_ste(1:length(ns)),'color',gray1,'linewidth',2)
errorbar(ns,keeps(3).Lhat_avg(1:length(ns)),keeps(3).Lhat_ste(1:length(ns)),'color',gray2,'linewidth',2)
axis([0 ns(end) 0 0.5])
xlabel('# training samples','fontsize',fs)
ylabel('misclassification rate','fontsize',fs)
% legend(alg(1).name,alg(2).name, alg(3).name,'fontsize',fs)
legend('coh','inc','nb'); %,'fontsize',fs,'Location','Best')
set(gca,'YTick',[0.1:0.1:0.5])



if savestuff==1
    print_fig(['../../figs/' fname '_Lhats'],[3 3]*1.5)
end

