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
        nstars=min(length(find(vscore>0)),num_stars);
        blank(vstars(1:nstars),:) = SigMat(vstars(1:nstars),:);
        blank(:,vstars(1:nstars))= SigMat(:,vstars(1:nstars));
        [~, indsp] = sort(blank(:));
        
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