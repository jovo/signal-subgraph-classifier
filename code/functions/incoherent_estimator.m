function ind = incoherent_estimator(SigMat,num_edges)
% finds the most significant num_edges edges  
% if num_edges isn't specified, it is set to # of vertices

if nargin==1, num_edges=sqrt(numel(SigMat)); end 
[~, delind]=sort(SigMat(:),'ascend');
ind  = delind(1:num_edges);