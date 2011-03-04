function [egg wcounter] = egg_estimator(SigMat,num_stars,num_edges)
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
wset=unique([0; sort(SigMat(:))]);
wset(wset>1-1e-3)=[];
wcounter = 1;
[V ~] = size(SigMat); 

noegg=0;
wconv=0;
inner_conv=0;
while wconv==0
    
    w=wset(wcounter);
    blank=SigMat;
    blank(blank>w)=1;
    blank(blank<=w)=0;
    score=V*2-(sum(blank)+sum(blank,2)');
    [vscore, vstars] = sort(score,'descend');
    
    if sum(vscore(1:num_stars))>=num_edges
        
        while inner_conv==0;
        
            blank=0*SigMat+1;
            nstars=min(length(find(vscore>0)),num_stars);
            blank(vstars(1:nstars),vstars(1:nstars)) = SigMat(vstars(1:nstars),vstars(1:nstars));

            blank(blank>w)=1;
            stilde=length(find((blank(:)<1)));
            
            if stilde>=num_edges
                inner_conv=1;
            else
                wcounter=wcounter+1;                
                w=wset(wcounter);
                blank=SigMat;
                blank(blank>w)=1;
                blank(blank<=w)=0;
                score=V*2-(sum(blank)+sum(blank,2)');
                [vscore, vstars] = sort(score,'descend');

                if wcounter==length(wset)
                    inner_conv=1;
                    noegg=1;
                end
            end
        
        end 
        
        [~, indsp] = sort(blank(:));
        
        egg=indsp(1:num_edges);
        wconv=1;
        
        if noegg==1, egg=NaN; end % if we never actually found an egg
    
    else
        wcounter=wcounter+1;
        if wcounter>length(wset),
            egg=[];
            wconv=1;
        end
    end % if statement
    
end % while state
wcounter=wcounter-1;
