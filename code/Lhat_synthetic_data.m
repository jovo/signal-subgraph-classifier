%% load real data
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


fname='synthetic_BLSA0317_FA';

%% generate synthetic data

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



%% alg stuff
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

%% classify
Out = plugin_classifier_cv_loop(As,ClassIDs,alg,'loo');


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
semilogx(Out(2).Lhat,'color','k','linewidth',2)
xlabel('log size of signal subgraph','fontsize',fs)
ylabel('misclassification rate','fontsize',fs)
title('incoherent estimator')
axis([1 1000 0 0.5])

h(2)=subplot(122);
imagesc(Out(1).Lhat), colorbar
xlabel('size of signal subgraph','fontsize',fs)
ylabel('number of star-vertices','fontsize',fs)
title('coherent estimator')

set(h,'fontsize',fs)
%%
print_fig(['../figs/' fname '_Lhats'],[8 3]*1.5)


