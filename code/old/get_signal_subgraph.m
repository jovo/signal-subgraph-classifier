clear, clc, clf
load('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/BLSA_0317_FAMtx');
constants = get_constants(AdjMats,ClassIDs);     % get constants to ease classification code
n=constants.s;

As=0*AdjMats;
for i=1:n
    A=AdjMats(:,:,i);
    A(A<0.2)=0;
    A(A>=0.2)=1;
    As(:,:,i)=A;
end


%%

V=constants.n;
num_star_list=1:V/2; star_len=length(num_star_list);
for m=1:star_len
    num_edge_list{m}=1:V/2*num_star_list(m);
end

phat    = get_ind_edge_params(As,constants);
SigMat  = run_get_fisher_pvals(As,constants);

incorrects=nan(star_len,length(num_edge_list{m}),n);

for m=1:star_len
    
    disp(num_star_list(m))
    edges=num_edge_list{m};
    
    for s=1:length(edges)
        
        for i=1:n % i indexes which graph we are leaving out
            coherent{m,s} = coherent_signal_subgraph_estimator(SigMat,num_star_list(m),edges(s));
            yhat = plugin_bern_classify(As,phat,coherent{m,s});
            incorrects(m,s,i)=yhat~=ClassIDs(i);
        end
    end
    save('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/get_signal_subgraph');
end

%% calc Lhat

siz=size(incorrects)
Lhat=nan(siz(1),siz(2));

for m=1:star_len
    for s=1:length(num_edge_list{m})    
        Lhat(m,s)=mean(incorrects(m,s,:));
    end
end


%% find optimal coherent signal subgraph

[foo ix] = min(Lhat(:));
[M S] = ind2sub(size(Lhat),ix);

blank=0*Sigmat;
blank(coherent{M,S})=1;

figure(1), clf
imagesc(blank)
title('Coherent Signal Subgraph')



%% find optimal star-1 signal subgraph

[foo ix] = min(Lhat(1,:));
[M S] = ind2sub(size(Lhat),ix);

blank=0*Sigmat;
blank(coherent{M,S})=1;

figure(2), clf
imagesc(blank)
title('Coherent Signal Subgraph')



save('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/get_signal_subgraph');
