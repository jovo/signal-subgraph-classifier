function [Lhat coherent incorrects] = coherent_classifier_loop(As,ys,star_list,edge_list,loo,Atst,ytst)
% loops over star_list and edge_list,
% computing the coherent signal subgraph and pluging classifier for each
% iteration.
%
% INPUTS:
%   As:         adjacency matrices
%   ys:         class labels
%   star_list:  list of how many star-vertices in coherent clasifier
%   edge_list:  array containing a list of the number of edges to use for each number of star-vertices
%   loo:        whether to do loo (optional: default is on)
%
% OUTPUTS:
%   Lhat:       average misclassification rate for each (m,s) pair
%   coherent:   coherent signal subgraph estimate
%   incorrects: whether the classifier got the answer correct

constants = get_constants(As,ys);
if nargin==4, loo=1; end

if loo==1
    
incorrects=nan(length(star_list),length(edge_list{end}),constants.s);
    for i=1:constants.s % i indexes which graph we are leaving out
        
        disp(['loocv iter: ' num2str(i)])
        
        [Atrn Gtrn Atst] = crossvalprep4(As,constants,[],i); % seperate data into training and testing sets
        phat    = get_ind_edge_params(Atrn,Gtrn);
        SigMat  = run_get_fisher_pvals(Atrn,Gtrn);
        
        for m=1:length(star_list)
            
            disp(['# star-vertics: ' num2str(star_list(m))])
            
            for s=1:length(edge_list{m})
                
                if mod(s,100)==0, disp(['# edges: ' num2str(edge_list{m}(s))]); end
                
                [coherent{m,s,i} incorrects(m,s,i)] = ...
                    coherent_classifier(SigMat, Atst, phat, star_list(m), edge_list{m}(s), ys(i));
                
            end
        end
    end
    
elseif loo==0 % learn signal subgraph from full data
    
incorrects=nan(length(star_list),length(edge_list{end}),constants.s);
    SigMat  = run_get_fisher_pvals(As,constants);
    
    for m=1:length(star_list)
        
        disp(['# star-vertics: ' num2str(star_list(m))])
        
        for s=1:length(edge_list{m})
            
            if mod(s,100)==0, disp(['# edges: ' num2str(edge_list{m}(s))]); end
            
            coherent{m,s} = coherent_estimator(SigMat,star_list(m),edge_list{m}(s));
            
            for i=1:constants.s % i indexes which graph we are leaving out
                [Atrn Gtrn Atst] = crossvalprep4(As,constants,[],i); % seperate data into training and testing sets
                phat    = get_ind_edge_params(Atrn,Gtrn);
                yhat = plugin_bern_classify(Atst,phat,coherent{m,s});
                incorrects(m,s,i)=yhat~=ys(i);
            end
        end
    end
    
elseif loo==2
    
    phat    = get_ind_edge_params(As,constants);
    SigMat  = run_get_fisher_pvals(As,constants);
    incorrects=nan(length(star_list),length(edge_list{end}),length(ytst));
    
    
    for i=1:length(ytst) % i indexes which graph we are leaving out
        
        disp(['test #: ' num2str(i)])
        
        for m=1:length(star_list)
            
            disp(['# star-vertics: ' num2str(star_list(m))])
            
            for s=1:length(edge_list{m})
                
                if mod(s,100)==0, disp(['# edges: ' num2str(edge_list{m}(s))]); end
                
                [coherent{m,s,i} incorrects(m,s,i)] = ...
                    coherent_classifier(SigMat, Atst, phat, star_list(m), edge_list{m}(s), ys(i));
                
            end
        end
    end
    
    
    
    
    
end


Lhat=nan(length(star_list),length(edge_list{end}));
for m=1:length(star_list)
    for s=1:length(edge_list{m})
        Lhat(m,s)=mean(incorrects(m,s,:));
    end
end