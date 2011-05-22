clear, clc, clf
load('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/BLSA_0317_FAMtx');

siz=size(AdjMats);
n=siz(3);           % # experiments
V=siz(1);           % # vertices

As=0*AdjMats;       % threshold
for i=1:n
    A=AdjMats(:,:,i);
    A(A<0.2)=0;
    A(A>=0.2)=1;
    As(:,:,i)=A;
end

star_list=1:V/2;    % list of # star-vertices
star_len=length(star_list);
for m=1:star_len    % list of # edges per star-vertex
    edge_list{m}=1:V*star_list(m)/2;
end

[Lhat coherent incorrects] = coherent_classifier_loop(As,ClassIDs,star_list,edge_list,1);

%% find best star_1 coherent signal subgraph using BLSA_0317 data

clc
star_list=1;
edge_list{1}=1:35;

[Lhat coherent incorrects] = coherent_classifier_loop(As,ClassIDs,star_list,edge_list,0);


%% find best star_1 coherent signal subgraph using BLSA50 data

clear, clc, clf
load('/Users/jovo/Research/data/MRI/BLSA/BLSA50/BLSA50_As_targs');

star_list=70;
edge_list{1}=1:1000;

[Lhat coherent incorrects] = coherent_classifier_loop(As,targs,star_list,edge_list,0);


%% synthetic
clear, clc

V=10;
m=1;
s=10;

p=0.5;
q=0.1;

Ess = [1 2 4 6 8 9];
thesem = zeros(V);
thesem(5,Ess) = 1;
thesem(Ess,5) = 1;
Ess=find(thesem);

E0=p*ones(V);
E1=p*ones(V);
E1(Ess)=q;

n=500;
As = nan(V,V,n);
As(:,:,1:n/2)=repmat(E0,[1 1 n/2]) > rand(V,V,n/2);
As(:,:,n/2+1:n)=repmat(E1,[1 1 n/2]) > rand(V,V,n/2);

ys=[zeros(1,n/2) ones(1,n/2)];


params.E0=E0;
params.E1=E1;

params.lnE0  = log(E0);
params.ln1E0 = log(1-E0);
params.lnE1  = log(E1);
params.ln1E1 = log(1-E1);

params.lnprior0 = log(1/2);
params.lnprior1 = log(1/2);
 
star_list=m;
edge_list{1}=s;


[Lhat coherent incorrects] = coherent_classifier_loop(As,ys,star_list,edge_list,1);

%% synthetic
clear, clc

V=10;
m=1;
s=10;

p=0.5;
q=0.1;

Ess = [1 2 4 6 8 9];
thesem = zeros(V);
thesem(5,Ess) = 1;
thesem(Ess,5) = 1;
Ess=find(thesem);

E0=p*ones(V);
E1=p*ones(V);
E1(Ess)=q;

n=500;
As = nan(V,V,n);
As(:,:,1:n/2)=repmat(E0,[1 1 n/2]) > rand(V,V,n/2);
As(:,:,n/2+1:n)=repmat(E1,[1 1 n/2]) > rand(V,V,n/2);

ys=[zeros(1,n/2) ones(1,n/2)];


params.E0=E0;
params.E1=E1;

params.lnE0  = log(E0);
params.ln1E0 = log(1-E0);
params.lnE1  = log(E1);
params.ln1E1 = log(1-E1);

params.lnprior0 = log(1/2);
params.lnprior1 = log(1/2);
 
alg.star_list=m;
alg.edge_list{1}=s;
cv='loo';


%%
[Lhat coherent incorrects] = plugin_classifier_cv_loop(As,ys,alg,cv);



%%


% params.d_pos = abs(E0-E1);           % position difference
% params.d_opt = abs(E0./sqrt(E0.*(1-E0)) - E1./sqrt(E1.*(1-E1))); % optimal difference
