function [Atrn Gtrn Atst Gtst inds] = crossvalprep4(As,constants,xval,i)
% this function prepares data for cross-validation using either:
% (1) hold-out data, or
% (2) using all permutations of leave-two-out (one from each class)
% depending on whether num_iters is sufficiently small to do all possible
% permuations
%
% INPUT:
%   As:         adjacency matrices
%   constants:  number of samples, etc.
%   xval:       structure containing xval parameters including
%       s0_trn: # of training samples for class 0
%       s1_trn: # of training samples for class 1
%       num_iters: (optional) max # of iters to perform
%       s0_tst: (optional) # of testing samples for class 0
%       s1_tst: (optional) # of testing samples for class 1
%
% OUTPUT:
%   Atrn:       adjacency matrices for training
%   Gtrn:       constatnts for training data
%   Atst:       adj. mat.'s for testing
%   Gtst:       constants for testing data
%   inds:       collection of indices (for graph invariant approach)


if ~isfield(xval,'loo')
    xval.loo=1;
end
if xval.loo==0;
    if ~isfield(xval,'s0_trn'), error('must specify xval.s0_trn'); end
    if ~isfield(xval,'s1_trn'), error('must specify xval.s1_trn'); end
    if ~isfield(xval,'num_iters'), error('must specify xval.num_iters'); end
end

siz=size(constants.y0); if siz(1)==1, y0=constants.y0; else y0=constants.y0'; end
siz=size(constants.y1); if siz(1)==1, y1=constants.y1; else y1=constants.y1'; end


%% select training and testing data

if xval.loo==true
    y0trn=y0;
    y1trn=y1;
    if any(i==y0)
        y0trn(y0==i)=[];
        y0tst=i;
        y1tst=[];
    else
        y1trn(y1==i)=[];
        y1tst=i;
        y0tst=[];
    end
else % randomly sample s0_trn and s1_trn data points for training, and use the others as testing
    if ~isfield(xval,'s0_tst'), xval.s0_tst=constants.s0-xval.s0_trn; end
    if ~isfield(xval,'s1_tst'), xval.s1_tst=constants.s1-xval.s1_trn; end
    
    ind0  = randperm(constants.s0);
    ind1  = randperm(constants.s1);
    
    y0trn = y0(ind0(1:xval.s0_trn));
    y0tst = y0(ind0(xval.s0_trn+1:xval.s0_trn+xval.s0_tst));
    
    y1trn = y1(ind1(1:xval.s1_trn));
    y1tst = y1(ind1(xval.s1_trn+1:xval.s1_trn+xval.s1_tst));
end


%% generate output

inds.ytrn=[y0trn y1trn];
inds.ytst=[y0tst y1tst];

Atrn = As(:,:,inds.ytrn);
Atst = As(:,:,inds.ytst);

ytrn = [zeros(1,length(y0trn)) ones(1,length(y1trn))];
ytst = [zeros(1,length(y0tst)) ones(1,length(y1tst))];

Gtrn = get_constants(Atrn,ytrn);
Gtst = get_constants(Atst,ytst);

if nargout==2
    inds.y0trn = y0trn;
    inds.y1trn = y1trn;
    inds.y0tst = y0tst;
    inds.y1tst = y1tst;
    inds.s_tst = length(inds.ytst);
end