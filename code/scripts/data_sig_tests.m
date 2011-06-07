%% load real data
clear, clc
load('~/Research/data/MRI/BLSA/BLSA_0317/base/BLSA_0317_countMtx');
fname='BLSA0317_Count_Lhats';
t=200;
siz=size(AdjMats);
n=siz(3);           % # experiments
V=siz(1);           % # vertices

As=0*AdjMats;       % threshold
for i=1:n
    A=AdjMats(:,:,i);
    A(A<t)=0;
    A(A>=t)=1;
    As(:,:,i)=tril(A,-1);
end


%% permute ClassIDs n_MC times

n_MC=10000;
parfor i=1:n_MC
    ix=randperm(n);
    Out = plugin_classifier_cv_loop(As,ClassIDs(ix),alg,'loo');
    Lhat(i)=Out.Lhat;   
end

save('../../data/data_sig_tests')


%% generate tests using estimated hyperparameters from cross-validation

i=0;
i=i+1;
alg(i).name='naive bayes';
alg(i).edge_list=find(tril(ones(V)-diag(ones(V,1)),-1));

i=i+1;
alg(i).name='incoherent';
alg(i).edge_list=10; %round(logspace(1,log10(nchoosek(V,2)/2),100));

i=i+1;
alg(i).name='coherent';
alg(i).star_list=12;
alg(i).edge_list{1}=360;

Out = plugin_classifier_cv_loop(As,ClassIDs,alg,'loo');


%% compare algs with chance
for i=1:3
    pval_chance=1-length(find(Out(i).Lhat<=Lhat))/n_MC;
end


%% mcnemar's test
clc

% nb vs. inc
b=double(~Out(1).incorrects & Out(2).incorrects);
c=double(Out(1).incorrects & ~Out(2).incorrects);

p_nb_inc=myBinomTest(sum(b),sum(b+c),0.5,'Two')

% nb vs. coh
b=double(~Out(1).incorrects & squeeze(Out(3).incorrects)');
c=double(Out(1).incorrects & ~squeeze(Out(3).incorrects)');

p_nb_coh=myBinomTest(sum(b),sum(b+c),0.5,'Two')

% inc vs. coh
b=double(~Out(2).incorrects & squeeze(Out(3).incorrects)');
c=double(Out(2).incorrects & ~squeeze(Out(3).incorrects)');

p_inc_coh=myBinomTest(sum(b),sum(b+c),0.5,'Two')


