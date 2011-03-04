%% synthetic
clear, clc

% constants
V=10;

% signal subgraph
m=1;
s=8;
Ess = 2:s+1;

% parameters
p=0.1;
q=0.2;

% parameters
E0=p*ones(V);
E1=p*ones(V);
E1(Ess)=q;

fname=['RE3_s' num2str(s) '_p' num2str(round(p*100)) '_q' num2str(round(q*100))];

%% simulate data for coherogram plots
clear h SigMat coherent incoherent coherogram Coherogram wset


ns=[10 50 100:100:500 1000];
nTrials=50;

inc_cor = nan(length(ns),1);
coh_cor = inc_cor;

for i=1:length(ns)
    
    disp(['nTrials = ' num2str(ns(i))])
    
    
    for t=1:nTrials
        
        % training data
        ntrn=ns(i);
        Atrn = nan(V,V,ntrn);
        Atrn(:,:,1:2:ntrn)=repmat(E0,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
        Atrn(:,:,2:2:ntrn)=repmat(E1,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
        for nn=1:ntrn, Atrn(:,:,nn)=tril(Atrn(:,:,nn),-1); end
        ytrn=repmat([0 1],1,ntrn/2);
        
        constants   = get_constants(Atrn,ytrn);     % get constants to ease classification code
        SigMat   = run_get_fisher_pvals(Atrn,constants);
        incoherent = incoherent_estimator(SigMat,s);
        [coherent cocount]= coherent_estimator(SigMat,m,s);
        
        inc_cor(i,t) = numcorrect(incoherent,Ess);
        coh_cor(i,t) = numcorrect(coherent,Ess);
        
    end
end

save(['../data/' fname])


%% get stats

% load(['../data/' fname])

inc_avg=mean(inc_cor,2);
coh_avg=mean(coh_cor,2);

inc_ste=std(inc_cor,[],2)./sqrt(nTrials);
coh_ste=std(coh_cor,[],2)./sqrt(nTrials);

eps=0.5;
rel_cor=inc_cor./(coh_cor+eps);

RE_avg=mean(rel_cor,2);
RE_ste=std(rel_cor')/sqrt(nTrials);
q=3;


figure(V), clf, clear h
grays=linspace(0,0.75,length(ns));
gray1=0.5*[1 1 1];
sa=0.5;
fs=10;

h(1)=subplot(121);
hold all
errorbar(ns,1-inc_avg/s,inc_ste/sqrt(nTrials),'-','linewidth',2,'color',gray1)
errorbar(ns,1-coh_avg/s,coh_ste/sqrt(nTrials),'-','linewidth',2,'color','k')
% title('edges','fontsize',fs)
ylabel('missed edge rate','fontsize',fs)
xlabel('# samples','fontsize',fs)
lh=legend('incoherent','coherent');
axis('tight')
set(lh,'Location','Best','fontsize',fs)
axis('tight')
set(gca,'fontsize',fs)

h(2)=subplot(122); hold all
errorbar(ns(q:end),RE_avg(q:end),RE_ste(q:end),'-','linewidth',2,'color','k')
plot([0 1000],[1 1],'--k')
ylabel('relative efficiency','fontsize',fs)
xlabel('# samples','fontsize',fs)
axis('tight')

set(h,'fontsize',fs)

%%
print_fig(['../figs/' fname '_RE'],[5 2]*2)

