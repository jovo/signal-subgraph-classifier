fname='synthetic_BLSA0317';

clear, clc
load('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/BLSA_0317_FAMtx');

siz=size(AdjMats);
n=siz(3);           % # experiments
V=siz(1);           % # vertices

t=0.4;
As=0*AdjMats;       % threshold
for i=1:n
    A=AdjMats(:,:,i);
    A(A<t)=0;
    A(A>=t)=1;
    As(:,:,i)=A;
end


% sufficient stats
constants   = get_constants(As,ClassIDs);
phat        = get_ind_edge_params(As,constants);
SigMat      = run_get_fisher_pvals(As,constants);

% parameters 
m=5;
s=20;
coherent = coherent_estimator(SigMat,m,s);

E0 = mean(As,3);
E1 = E0;
E0(coherent)=phat.E0(coherent);
E1(coherent)=phat.E1(coherent);

% synethtic graphs
As = nan(V,V,n);
As(:,:,1:constants.s0)=repmat(E0,[1 1 constants.s0]) > rand(V,V,constants.s0);
As(:,:,constants.s0+1:end)=repmat(E1,[1 1 constants.s1]) > rand(V,V,constants.s1);

ytrn=[zeros(constants.s0,1); ones(constants.s1,1)];

fname='synthetic_BLSA0317_rates';

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
ns=[10:10:50 100:50:300];


for i=1:length(ns)
    
    disp(['ntrn ', num2str(ns(i))])

    for j=1:nTrials
        
        disp(['trial ', num2str(j)]) 
        % training data
        ntrn=ns(i);
        Atrn = nan(V,V,ntrn);
        Atrn(:,:,1:2:ntrn)=repmat(E0,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
        Atrn(:,:,2:2:ntrn)=repmat(E1,[1 1 ntrn/2]) > rand(V,V,ntrn/2);

        ytrn=repmat([0 1],1,ntrn/2);
        
        % testing data
        ntst=100;
        Atst = nan(V,V,ntst);
        Atst(:,:,1:2:ntst)=repmat(E0,[1 1 ntst/2]) > rand(V,V,ntst/2);
        Atst(:,:,2:2:ntst)=repmat(E1,[1 1 ntst/2]) > rand(V,V,ntst/2);
        ytst=repmat([0 1],1,ntst/2);
        
        % classify
        Outtemp = plugin_classifier_cv_loop(Atrn,ytrn,alg,'HoldOut',Atst,ytst);
        
        % get stats
        Out(1).Lhat(i,j)   = Outtemp(1).Lhat;
        Out(1).correct_edges(i,j) = numcorrect(Outtemp(1).coherent{1},coherent);
        
        Out(2).Lhat(i,j)   = Outtemp(2).Lhat;
        Out(2).correct_edges(i,j) = numcorrect(Outtemp(2).incoherent{1},coherent);
        
        Out(3).Lhat(i,j)   = Outtemp(3).Lhat;
    end
    
end
save(['../data/' fname])

%% plot rates

load(['../data/' fname])

figure(1), clf

gray1=0.5*[1 1 1];
gray2=0.75*[1 1 1];
sa=0.5;
fs=10;

h(1)=subplot(121); hold all
errorbar(ns,mean(Out(1).correct_edges,2),std(Out(1).correct_edges,[],2)/sqrt(nTrials),'k','linewidth',2)
errorbar(ns,mean(Out(2).correct_edges,2),std(Out(2).correct_edges,[],2)/sqrt(nTrials),'color',gray1,'linewidth',2)
axis([0 ns(end) 0 s+0.1])
set(gca,'YTick',[2:2:s])
xlabel('# training samples','fontsize',fs)
ylabel('# identified edges','fontsize',fs)

h(2)=subplot(122); hold all
errorbar(ns,mean(Out(1).Lhat,2),std(Out(1).Lhat,[],2)/sqrt(nTrials),'k','linewidth',2)
errorbar(ns,mean(Out(2).Lhat,2),std(Out(2).Lhat,[],2)/sqrt(nTrials),'color',gray1,'linewidth',2)
errorbar(ns,mean(Out(3).Lhat,2),std(Out(3).Lhat,[],2)/sqrt(nTrials),'color',gray2,'linewidth',2)
axis([0 ns(end) 0 0.5])
xlabel('# training samples','fontsize',fs)
ylabel('misclassification rate','fontsize',fs)
% legend(alg(1).name,alg(2).name, alg(3).name,'fontsize',fs)
legend('coh','inc','nb','fontsize',fs,'Location','Best')
set(gca,'YTick',[0.1:0.1:0.5])



for i=1:numel(h)
    op=get(h(i),'OuterPosition');
    op(4)=op(4)+op(2);
    op(2)=0;
    set(h(i),'OuterPosition',op)
    set(h(i),'fontsize',fs)
end

print_fig(['../figs/' fname '_rates'],[6 2]*1.5)


