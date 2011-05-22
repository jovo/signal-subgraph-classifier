%% generate parameters for synthetic homo data

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

star_list=1%:V/2;    % list of # star-vertices
star_len=length(star_list);
for m=1:star_len    % list of # edges per star-vertex
    edge_list{m}=1:V*star_list(m)/1;
end

[Lhat coherent incorrects] = coherent_classifier_loop(As,ClassIDs,star_list,edge_list,0);

%%

[foo min_id] = min(Lhat);

SigSubGraph = zeros(V);
SigSubGraph(coherent{1,min_id})=1;
InSigInds=find(triu(SigSubGraph==0,1));
InSigSubGraph = zeros(V);
InSigSubGraph(InSigInds)=1;


SigVec=nan(min_id,n);
InSigVec=nan(numel(InSigInds),n);
for i=1:n
    A=As(:,:,i);
    SigVec(:,i)=A(coherent{1,min_id});
    InSigVec(:,i)=A(InSigInds);    
end

p=mean(SigVec(:));
q=mean(InSigVec(:));

subplot(121), imagesc(SigVec), title(p)
subplot(122), imagesc(InSigVec), title(q)


InSig1=nan(numel(InSigInds),length(find(ClassIDs)));
for i=find(ClassIDs)
    A=As(:,:,i);
    InSig1(:,i)=A(InSigInds);    
end
InSig1=InSig1(:,find(ClassIDs));
q1=mean(InSig1(:));

InSig0=nan(numel(InSigInds),length(find(~ClassIDs)));
for i=find(~ClassIDs)
    A=As(:,:,i);
    InSig0(:,i)=A(InSigInds);    
end
InSig0=InSig0(:,find(~ClassIDs));
q0=mean(InSig0(:));


%%

M.n=V;
M.s=1000;

M.E0=p*tril(ones(V),-1);
M.E1=p*ones(V);
M.E1(coherent{1,min_id})=q;
M.E1=tril(M.E1',-1);

sim = gen_ind_edge_data(M)

save('/Users/jovo/Research/data/MRI/BLSA/synthetic/homo.mat','M','sim','As','ClassIDs');


%%

p=0.5;
qd=0.2;
q0=p-qd;
q1=p+qd;

V=10;
M.n=V;
M.s=1000;

SigID=2:8;

M.E0=p*tril(ones(V),-1);
M.E0(SigID)=q0;
M.E0=tril(M.E0,-1);

M.E1=p*ones(V);
M.E1(SigID)=q1;
M.E1=tril(M.E1,-1);

sim = gen_ind_edge_data(M)

save('/Users/jovo/Research/data/MRI/BLSA/synthetic/homo2.mat','M','sim');


%%

p=0.1;
q=0.5;

V=279;
M.n=V;
M.s=520;
m=36;

thism=1;
thesem=2:m+1;

M.E0=p*tril(ones(V),-1);

M.E1=p*ones(V);
M.E1(thesem,thism)=q;
M.E1=tril(M.E1,-1);

sim = gen_ind_edge_data(M)

save('/Users/jovo/Research/data/MRI/BLSA/synthetic/celegans.mat','M','sim');



%% homo 1

p=0.1;
q=0.5;

V=10;
M.n=V;
M.s=500;
m=5;

thism=1;
thesem=2:m+1;

M.E0=p*tril(ones(V),-1);

M.E1=p*ones(V);
M.E1(thesem,thism)=q;
M.E1=tril(M.E1,-1);

trn = gen_ind_edge_data(M)

M.s=20;
tst = gen_ind_edge_data(M)

save('/Users/jovo/Research/data/MRI/BLSA/synthetic/homo1.mat','trn','tst','M');
