function [egg wcounter] = egg2_estimator(SigMat,num_stars,num_edges)
% estimates the coherent signal subgraph from
% INPUTS:
%   SigMat:     matrix of significance of differences
%   num_stars:  # of star-vertices
%   num_edges:  # of edges
%
% OUTPUT:
%   egg:   signal subgraph estimate
%   wcounter:   how many levels we had to extend to get enough edges

%% code
wset=unique(sort(SigMat(:)));
wset(wset>1-1e-3)=[];
wcounter = 1;
siz = size(SigMat);
V = siz(1);

wconv=0;
while wconv==0
    
    w=wset(wcounter);
    inds = find(SigMat<=w);
    [I J] = ind2sub(siz,inds);
    ncounts = histc([I; J],1:V);
    [B IX] = sort(ncounts,'descend');
    sumcounts = sum(B(1:num_stars));
    if sumcounts>=num_edges
        blank=0*SigMat+1;
        blank(IX(1:num_stars),IX(1:num_stars))=SigMat(IX(1:num_stars),IX(1:num_stars));
        [~, indsp] = sort(blank(:));
        egg=indsp(1:num_edges);
        wconv=1;
    else
        wcounter=wcounter+1;
        if wcounter>numel(wset),
            blank=0*SigMat+1;
            blank(IX(1:num_stars),IX(1:num_stars))=SigMat(IX(1:num_stars),IX(1:num_stars));
            [~, indsp] = sort(blank(:));
            egg=indsp(1:num_edges);
            wconv=1;
        end
    end
    
end % while wconv
