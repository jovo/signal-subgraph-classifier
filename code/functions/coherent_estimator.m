function [coherent wset wcounter] = coherent_estimator(SigMat,num_stars,num_edges)
% estimates the coherent signal subgraph from
% INPUTS:
%   SigMat:     matrix of significance of differences
%   num_stars:  # of star-vertices
%   num_edges:  # of edges
%
% OUTPUT:
%   coherent:   signal subgraph estimate
%   wcounter:   how many levels we had to extend to get enough edges

%% code
wset=unique([0; sort(SigMat(:))]);
% wset(wset>1-1e-3)=[];
wcounter = 1;
[V ~] = size(SigMat); 

wconv=0;
while wconv==0
    
    w=wset(wcounter);
    blank=SigMat;
    blank(blank>w)=1;
    blank(blank<=w)=0;
    score=V*2-(sum(blank)+sum(blank,2)');
    [vscore, vstars] = sort(score,'descend');
    
    if sum(vscore(1:num_stars))>=num_edges
        
        blank=0*SigMat+1;
        num_candidate_stars=length(find(vscore>0));
        if num_candidate_stars>num_stars
            worst_candidates=find(vscore==vscore(num_stars));
            best_candidates=find(vscore<vscore(num_stars));
            worst_ind=randperm(length(worst_candidates));
            worst_n=length(num_stars)-length(best_candidates);
            vhat_stars=[best_candidates worst_candidates(worst_ind(1:worst_n))];            
        else 
            vhat_stars=vstars(1:nstars);
        end
        blank(vhat_stars,:) = SigMat(vhat_stars,:);
        blank(:,vhat_stars)= SigMat(:,vhat_stars);
        [foo, indsp] = sort(blank(:));
        
        % keep all edges that are as significant as those that are counting
        %         coherent=indsp(1:find(foo>foo(num_edges),1));
        coherent=indsp(1:num_edges);
        
        wconv=1;
    else
        wcounter=wcounter+1;
        if wcounter>length(wset),
            coherent=[];
            wconv=1;
        end
    end % if statement
    
end % while state
wcounter=wcounter-1;