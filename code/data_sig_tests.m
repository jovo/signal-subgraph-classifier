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

%%
mark=[
0	0	1	0.062244753606925635
1	1	1	0.5483292201040381
2	1	1	1.2234201546825252
3	0	0	0.040089061874994766
4	1	1	0.45417776126056264
5	0	0	0.2628290021243338
6	0	0	0.2909978992578783
7	0	0	0.3132130786385885
8	1	0	0.3309547488756656
9	0	0	0.11041403480804415
10	0	0	0.05383168044193106
11	0	0	0.07129005543551464
12	1	1	0.24279790655592232
13	1	1	0.15020905992522357
14	0	0	0.03297778335180895
15	1	1	0.1023409968010022
16	0	0	0.1924228194096547
17	0	0	0.04264939303611927
18	0	0	0.06088115938687777
19	0	0	0.36582882940686984
20	1	1	0.3271651053215823
21	0	0	0.017655818204648753
22	1	1	0.07155107265150827
23	0	0	0.11895104180850076
24	1	0	0.04565490392644502
25	1	0	0.18294223601840434
26	0	1	0.1626189464522493
27	0	1	0.09871863119568819
28	0	1	0.10426584010368209
29	1	1	0.17859216981696896
30	1	1	0.04218313392038722
31	1	0	0.15721071981059354
32	0	1	0.031259635398005206
33	1	1	0.7656458772919794
34	1	1	0.0523223038682149
35	1	1	0.6560725428193254
36	0	0	0.031002458860815502
37	0	0	0.06695083347470061
38	0	0	0.002453961407115841
39	1	1	0.06247476896743547
40	0	1	0.035321294963294975
41	1	1	0.21401653463910048
42	0	1	0.5585967496007
43	0	1	0.060032062309077705
44	1	0	0.018630510959995978
45	0	0	0.20524636426234616
46	1	1	0.1928348215889144
47	0	0	0.04226982086126173
48	1	1	0.0610230585397991];

mark_wrong=mark(:,2)~=mark(:,3);
Out(4).Lhat=sum(mark_wrong)/length(mark_wrong);


%% compare algs with chance
for i=1:4
    pval_chance(i)=1-length(find(Out(i).Lhat<=Lhat))/n_MC;
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


save('../../data/data_sig_tests')

% mark vs. coh
b=double(~mark_wrong & squeeze(Out(3).incorrects));
c=double(mark_wrong & ~squeeze(Out(3).incorrects));

p_mark_coh=myBinomTest(sum(b),sum(b+c),0.5,'Two')

% mark vs. inc
b=double(~mark_wrong & (Out(2).incorrects)');
c=double(mark_wrong & ~(Out(2).incorrects)');

p_mark_inc=myBinomTest(sum(b),sum(b+c),0.5,'Two')

% mark vs. nb
b=double(~mark_wrong & (Out(1).incorrects)');
c=double(mark_wrong & ~(Out(1).incorrects)');

p_mark_nb=myBinomTest(sum(b),sum(b+c),0.5,'Two')


