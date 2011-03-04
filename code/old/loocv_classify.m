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
num_star_list=1:V/2;
star_len=length(num_star_list);
for m=1:star_len
    num_edge_list{m}=1:V/2*num_star_list(m);
end

incorrects=nan(star_len,length(num_edge_list{m}),n);


for i=1:n % i indexes which graph we are leaving out
    
    disp(['loocv iter: ' num2str(i)])
    
    [Atrn Gtrn Atst Gtst inds] = crossvalprep4(As,constants,[],i); % seperate data into training and testing sets
    phat    = get_ind_edge_params(Atrn,Gtrn);
    SigMat  = run_get_fisher_pvals(Atrn,Gtrn);
    
    for m=1:star_len
        
        disp(num_star_list(m))
        edges=num_edge_list{m};
        
        for s=1:length(edges)
            
            temp = coherent_signal_subgraph_estimator(SigMat,num_star_list(m),edges(s));
            
            if ~isempty(temp)
                coherent{m,s}=temp;
                yhat = plugin_bern_classify(Atst,phat,coherent{m,s});
                incorrects(m,s,i)=yhat~=ClassIDs(i);
            else
                break
            end
            
        end
    end
    save('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/loocv_classify');
end

%% calc Lhat

Lhat=nan(star_len,length(num_edge_list));

for m=1:star_len
    for s=1:length(num_edge_list{m})
        Lhat(m,s)=mean(incorrects(m,s,:));
    end
end

save('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/loocv_classify');
