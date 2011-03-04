function [coherent incorrect] = coherent_classifier(SigMat, AdjMat, phat, num_stars, num_edges, ClassID)
% implements the coherent classifier
% INPUTS:
%   SigMat:     matrix of significances
%   num_stars:  # of star-vertices
%   num_edges:  # of edges
%   ClassID:    true class (when available)
% OUTPUTS:
%   coherent:   the coherent signal subgraph
%   incorrect:  whether yhat=y


coh = coherent_estimator(SigMat,num_stars,num_edges);

if ~isempty(coh)
    coherent=coh;
    yhat = plugin_bern_classify(AdjMat,phat,coherent);
    if nargin==6, incorrect=yhat~=ClassID; else incorrect=NaN; end
else
    coherent=[];
    incorrect=NaN;
end