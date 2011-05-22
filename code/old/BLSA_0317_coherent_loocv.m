clear, clc, clf
load('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/BLSA_0317_FAMtx');

siz=size(AdjMats);
n=siz(3);           % # experiments
V=siz(1);           % # vertices

star_list=1%:V/2;    % list of # star-vertices
star_len=length(star_list);
for m=1:star_len    % list of # edges per star-vertex
    edge_list{m}=1:V/2*star_list(m);
end

alg(1).name='incoherent';
alg(1).edge_list=unique(round(logspace(0,log10(2415),100)));

ts=[0.2 0.3 0.4 0.5];

for t=1:length(ts);
    As=0*AdjMats;       % threshold
    for i=1:n
        A=AdjMats(:,:,i);
        A(A<ts(t))=0;
        A(A>=ts(t))=1;
        As(:,:,i)=A;
    end
    
    Out1 = plugin_classifier_cv_loop(As,ClassIDs,alg,'loo');
    Lhat{t}=Out1.Lhat;
    
end

%%
figure(1), clf, hold all
for i=1:length(ts)
    
    semilogx(Lhat{i})
end

%%

clear alg
alg(1).name='coherent';
alg(1).star_list=V/2:1:V;
for m=1:length(alg(1).star_list)    % list of # edges per star-vertex
    alg(1).edge_list{m}=unique(round(logspace(0,log10(V/2*alg(1).star_list(m)),100)));
end

%%
ts=0.4; %[0.2 0.3 0.4 0.5];

As=0*AdjMats;       % threshold
for i=1:n
    A=AdjMats(:,:,i);
    A(A<ts)=0;
    A(A>=ts)=1;
    As(:,:,i)=A;
end

Out2 = plugin_classifier_cv_loop(As,ClassIDs,alg,'loo');
Lhat{t}=Out.Lhat;

%%

constants = get_constants(As,ClassIDs);
SigMat  = run_get_fisher_pvals(As,constants);

%%

clear alg
alg(1).name='incoherent';
alg(1).edge_list=1:5;

ts=0.4; %[0.2 0.3 0.4 0.5];

As=0*AdjMats;       % threshold
for i=1:n
    A=AdjMats(:,:,i);
    A(A<ts)=0;
    A(A>=ts)=1;
    As(:,:,i)=A;
end

Out3 = plugin_classifier_cv_loop(As,ClassIDs,alg,'loo');

%%

clear alg
alg(1).name='coherent';
alg(1).star_list=35;
alg(1).edge_list{1}=1:5;

ts=0.4; %[0.2 0.3 0.4 0.5];

As=0*AdjMats;       % threshold
for i=1:n
    A=AdjMats(:,:,i);
    A(A<ts)=0;
    A(A>=ts)=1;
    As(:,:,i)=A;
end

Out3 = plugin_classifier_cv_loop(As,ClassIDs,alg,'loo');

%%

clear alg
alg(1).name='coherent';
alg(1).star_list=1:V;
for m=1:length(alg(1).star_list)    % list of # edges per star-vertex
    alg(1).edge_list{m}=unique(round(logspace(0,log10(V/2*alg(1).star_list(m)),100)));
end

ts=0.4; %[0.2 0.3 0.4 0.5];

As=0*AdjMats;       % threshold
for i=1:n
    A=AdjMats(:,:,i);
    A(A<ts)=0;
    A(A>=ts)=1;
    As(:,:,i)=A;
end

Out4 = plugin_classifier_cv_loop(As,ClassIDs,alg,'loo');
Lhat{t}=Out.Lhat;

