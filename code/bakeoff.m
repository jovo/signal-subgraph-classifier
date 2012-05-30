%% load real data
clear all
load('../../data/base/BLSA0317');
fname='BLSA0317_bakeoff';

As=Abinary;
Ys=ClassIDs;
[V , ~, n] = size(As);



%% elastic net classify
[B,FitInfo] = lassoglm([ones(1, n);  reshape(As,[V^2 n])]',Ys','binomial','CV',10);
lassoPlot(B,FitInfo,'plottype','CV');

yhat=nan(n,1);
parfor k=1:n
    yhat(k) = glmval(B(:,FitInfo.IndexMinDeviance),reshape(As(:,:,k),[V^2 1])','logit');
end
Lhats.lasso=sum(yhat>0.75~=Ys')/n;
yhats.lasso = yhat>0.75~=Ys';
% best performing lasso: 0.2653

%% kNN classify
kvec=1:48;
labeledIDM = InterpointDistanceMatrix(As);
[Lhats.knn yhats.knn] = knnclassifyIDM(labeledIDM,Ys',kvec);
figure(2)
plot(kvec,Lhats.knn);
[minval minidxk] = min(Lhats.knn)
% minval = 0.3469


%% invariants classify

x = get_graph_invariants(As,[1:2, 4:7]);
x0=x-repmat(mean(x'),n,1)';
x1=x0./repmat(std(x0'),n,1)';
labeledIDM = InterpointDistanceMatrix(x1);
[Lhats.invariants yhats.invariants] = knnclassifyIDM(labeledIDM,Ys',kvec);
figure(3)
plot(kvec,Lhats.invariants)

[minval minidxi] = min(Lhats.invariants)
% minval = 0.4286
% minidxi = 12
    
%% save stuff

savestuff=1;
if savestuff==1, save(['../../data/' fname]); end

%% stats

load('../../data/results/data_sig_tests')

% nb vs. coh
b=double(~Out(1).incorrects & squeeze(Out(3).incorrects)');
c=double(Out(1).incorrects & ~squeeze(Out(3).incorrects)');
p_nb_coh=myBinomTest(sum(b),sum(b+c),0.5,'Lesser')

% inc vs. coh
b=double(~Out(2).incorrects & squeeze(Out(3).incorrects)');
c=double(Out(2).incorrects & ~squeeze(Out(3).incorrects)');
p_inc_coh=myBinomTest(sum(b),sum(b+c),0.5,'Lesser')

% lasso vs. coh
b=double(~(yhats.lasso~=Ys') & squeeze(Out(3).incorrects));
c=double(yhats.lasso~=Ys' & ~squeeze(Out(3).incorrects));
p_lasso_coh=myBinomTest(sum(b),sum(b+c),0.5,'Lesser')

% invariants vs. coh
b=double(~(yhats.invariants(minidxi,:)'~=Ys') & squeeze(Out(3).incorrects));
c=double(yhats.invariants(minidxi,:)'~=Ys' & ~squeeze(Out(3).incorrects));
p_inv_coh=myBinomTest(sum(b),sum(b+c),0.5,'Lesser')

% G-kNN vs. coh
b=double(~(yhats.knn(minidxk,:)'~=Ys') & squeeze(Out(3).incorrects));
c=double(yhats.knn(minidxk,:)'~=Ys' & ~squeeze(Out(3).incorrects));
p_knn_coh=myBinomTest(sum(b),sum(b+c),0.5,'Lesser')

