function [Out] = plugin_classifier_cv_loop(As,ys,alg,cv,Atst,ytst)
% loops over star_list and edge_list,
% computing the coherent signal subgraph and pluging classifier for each
% iteration.
%
% INPUTS:
%   As:         adjacency matrices
%   ys:         class labels
%   cv:         cross-validation type: 'loov','InSample','HoldOut'
%   alg:        structure containing algorithm parameters
%   Atst:       (only when cv='HoldOut'): test graphs
%   ytst:       (only when cv='HoldOut'): test class
%
% OUTPUTS:
%   Out:        structure with the following fields (one per algorithm)
%   Lhat:       average misclassification rates
%   incorrects: whether the classifier got the answer correct
%   other fields specific to each algoritm

constants = get_constants(As,ys);
if nargin==3, cv=1; end

%% set up stuff for different algorithms
coh=0; inc=0; nb=0;
for a=1:numel(alg)
    if strcmp(alg(a).name,'coherent'),
        coh=1;
        coh_id=a;
        
        if nargin>4, n_inc=length(ytst); else n_inc=length(ys); end
        Out(a).incorrects=nan(length(alg(a).star_list),length(alg(a).edge_list{end}),n_inc);
        
    elseif strcmp(alg(a).name,'egg'),
        egg=1;
        egg_id=a;
        
        if nargin>4, n_inc=length(ytst); else n_inc=length(ys); end
        Out(a).incorrects=nan(length(alg(a).star_list),length(alg(a).edge_list{end}),n_inc);

    elseif strcmp(alg(a).name,'incoherent'),
        inc=1;
        inc_id=a;
        
        if nargin>4, n_inc=length(ytst); else n_inc=length(ys); end
        Out(a).incorrects=nan(length(alg(a).edge_list),n_inc);
        
    elseif strcmp(alg(a).name,'naive bayes'),
        nb=1;
        nb_id=a;
    end
    Out(a).name=alg(a).name;
end


%% cross-validate over different algorithms
if strcmp(cv,'loo')    % get signal subgraph using training data, then classify test data
    
    
    for i=1:constants.s, disp(['loocv iter: ' num2str(i)])
        
        [Atrn Gtrn Atst] = crossvalprep4(As,constants,[],i); % seperate data into training and testing sets
        phat    = get_ind_edge_params(Atrn,Gtrn);
        
        if coh==1 || inc==1 || egg==1
            SigMat  = run_get_fisher_pvals(Atrn,Gtrn);
        end
        
        if coh==1 % iterate over m & s
            a=coh_id;
            for m=1:length(alg(a).star_list), disp(['# star-vertics: ' num2str(alg(a).star_list(m))])
                for s=1:length(alg(a).edge_list{m}), if mod(s,100)==0, disp(['# edges: ' num2str(alg(a).edge_list{m}(s))]); end
                    Out(a).subspace{m,s,i} = coherent_estimator(SigMat,alg(a).star_list(m),alg(a).edge_list{m}(s));
                    Out(a).incorrects(m,s,i) = plugin_bern_classify(Atst,phat,Out(a).subspace{m,s,i},ys(i));
                end
            end
        end
        
        if egg==1 % iterate over m & s
            a=egg_id;
            for m=1:length(alg(a).star_list), disp(['# star-vertics: ' num2str(alg(a).star_list(m))])
                for s=1:length(alg(a).edge_list{m}), if mod(s,100)==0, disp(['# edges: ' num2str(alg(a).edge_list{m}(s))]); end
                    Out(a).subspace{m,s,i} = egg2_estimator(SigMat,alg(a).star_list(m),alg(a).edge_list{m}(s));
                    Out(a).incorrects(m,s,i) = plugin_bern_classify(Atst,phat,Out(a).subspace{m,s,i},ys(i));
                end
            end
        end
        
        if inc==1 % iterate over s
            a=inc_id;
            edge_list=alg(inc_id).edge_list;
            for s=1:length(edge_list), if mod(s,100)==0, disp(['# edges: ' num2str(edge_list(s))]); end
                Out(inc_id).subspace{s,i} = incoherent_estimator(SigMat,edge_list(s));
                Out(inc_id).incorrects(s,i) = plugin_bern_classify(Atst,phat,Out(inc_id).subspace{s,i},ys(i));
            end
        end
        
        if nb==1 % naive bayes
            Out(nb_id).incorrects(i) = plugin_bern_classify(Atst,phat,alg(nb_id).edge_list,ys(i));
        end            
        
    end
    
    
elseif strcmp(cv,'InSample') % learn signal subgraph from full data
    
    if coh==1 || inc==1 || egg==1
        SigMat  = run_get_fisher_pvals(As,constants);
    end
    
    if coh==1 % iterate over m & s
        a=coh_id;
        for m=1:length(alg(a).star_list), 
            if numel(alg(a).star_list)>1, disp(['# star-vertics: ' num2str(alg(a).star_list(m))]); end
            for s=1:length(alg(a).edge_list{m}), if mod(s,100)==0, disp(['# edges: ' num2str(alg(a).edge_list{m}(s))]); end
                
                Out(a).subspace{m,s} = coherent_estimator(SigMat,alg(a).star_list(m),alg(a).edge_list{m}(s));
                for i=1:constants.s % i indexes which graph we are leaving out
                    [Atrn Gtrn Atst] = crossvalprep4(As,constants,[],i); % seperate data into training and testing sets
                    phat    = get_ind_edge_params(Atrn,Gtrn);
                    Out(a).incorrects(m,s,i) = plugin_bern_classify(Atst,phat,Out(a).subspace{m,s},ys(i));
                end
            end
        end
        
    end
    
    if egg==1 % iterate over m & s
        a=egg_id;
        for m=1:length(alg(a).star_list), 
            if numel(alg(a).star_list)>1, disp(['# star-vertics: ' num2str(alg(a).star_list(m))]); end
            for s=1:length(alg(a).edge_list{m}), if mod(s,100)==0, disp(['# edges: ' num2str(alg(a).edge_list{m}(s))]); end
                
                Out(a).subspace{m,s} = egg2_estimator(SigMat,alg(a).star_list(m),alg(a).edge_list{m}(s));
                for i=1:constants.s % i indexes which graph we are leaving out
                    [Atrn Gtrn Atst] = crossvalprep4(As,constants,[],i); % seperate data into training and testing sets
                    phat    = get_ind_edge_params(Atrn,Gtrn);
                    Out(a).incorrects(m,s,i) = plugin_bern_classify(Atst,phat,Out(a).subspace{m,s},ys(i));
                end
            end
        end
        
    end
    
    
    if inc==1 % iterate over s
        a=inc_id;
        for s=1:length(alg(inc_id).edge_list), if mod(s,100)==0, disp(['# edges: ' num2str(alg(inc_id).edge_list(s))]); end
            
            Out(inc_id).subspace{s} = incoherent_estimator(SigMat,alg(inc_id).edge_list(s));
            for i=1:constants.s % i indexes which graph we are leaving out
                [Atrn Gtrn Atst] = crossvalprep4(As,constants,[],i); % seperate data into training and testing sets
                phat    = get_ind_edge_params(Atrn,Gtrn);
                Out(inc_id).incorrects(s,i) = plugin_bern_classify(Atst,phat,Out(inc_id).subspace{s},ys(i));
            end
        end
        
    end
    
    if nb==1
        for i=1:constants.s % i indexes which graph we are leaving out
            [Atrn Gtrn Atst] = crossvalprep4(As,constants,[],i); % seperate data into training and testing sets
            phat    = get_ind_edge_params(Atrn,Gtrn);
            Out(nb_id).incorrects(i) = plugin_bern_classify(Atst,phat,alg(nb_id).edge_list,ys(i));
        end
    end
    
    
elseif strcmp(cv,'HoldOut'), % get signal subgraph from training data, test on held-out data
    
    phat    = get_ind_edge_params(As,constants);
    
    if coh==1 || inc==1 || egg==1
        SigMat  = run_get_fisher_pvals(As,constants);
    end
    
    if coh==1 % iterate over m & s
        a=coh_id;
        for m=1:length(alg(a).star_list), 
            if numel(alg(a).star_list)>1, disp(['# star-vertics: ' num2str(alg(a).star_list(m))]); end
            for s=1:length(alg(a).edge_list{m}), if mod(s,100)==0, disp(['# edges: ' num2str(alg(a).edge_list{m}(s))]); end
                
                Out(a).subspace{m,s} = coherent_estimator(SigMat,alg(a).star_list(m),alg(a).edge_list{m}(s));
                for i=1:length(ytst) % i indexes which test graph
                    Out(a).incorrects(m,s,i) = plugin_bern_classify(Atst(:,:,i),phat,Out(a).subspace{m,s},ytst(i));
                end
                
            end % s
        end % m
    end % coh

        
    if egg==1 % iterate over m & s
        a=egg_id;
        for m=1:length(alg(a).star_list), 
            if numel(alg(a).star_list)>1, disp(['# star-vertics: ' num2str(alg(a).star_list(m))]); end
            for s=1:length(alg(a).edge_list{m}), if mod(s,100)==0, disp(['# edges: ' num2str(alg(a).edge_list{m}(s))]); end
                
                Out(a).subspace{m,s} = egg2_estimator(SigMat,alg(a).star_list(m),alg(a).edge_list{m}(s));
                for i=1:length(ytst) % i indexes which test graph
                    Out(a).incorrects(m,s,i) = plugin_bern_classify(Atst(:,:,i),phat,Out(a).subspace{m,s},ytst(i));
                end
                
            end % s
        end % m
    end % egg

    
    if inc==1 % iterate over s
        a=inc_id;
        for s=1:length(alg(inc_id).edge_list), if mod(s,100)==0, disp(['# edges: ' num2str(alg(inc_id).edge_list(s))]); end
            
            Out(inc_id).subspace{s} = incoherent_estimator(SigMat,alg(inc_id).edge_list(s));
            for i=1:length(ytst) % i indexes which test graph
                Out(inc_id).incorrects(s,i) = plugin_bern_classify(Atst(:,:,i),phat,Out(inc_id).subspace{s},ytst(i));
            end % i
            
        end % s
    end % inc
    
    if nb==1 % naive bayes
        for i=1:length(ytst) % i indexes which test graph
            Out(nb_id).incorrects(i) = plugin_bern_classify(Atst(:,:,i),phat,alg(nb_id).edge_list,ytst(i));
        end % i
    end % nb
    
    
end % cv



%% compute Lhat for each algorithm
if coh==1
    Out(coh_id).Lhat=mean(Out(coh_id).incorrects,3);
end

if egg==1
    Out(egg_id).Lhat=mean(Out(egg_id).incorrects,3);
end

if inc==1
    Out(inc_id).Lhat=mean(Out(inc_id).incorrects,2);
end

if nb==1
    Out(nb_id).Lhat=mean(Out(nb_id).incorrects);
end