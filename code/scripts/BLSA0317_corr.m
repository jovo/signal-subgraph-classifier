clear, clc
load('../../data/BLSA0317_Count_Lhats');

siz=size(AdjMats);
n=siz(3);           % # experiments
V=siz(1);           % # vertices

t=0.4;
As=0*AdjMats;       % threshold
for i=1:n
    A=AdjMats(:,:,i);
    A(A<t)=0;
    A(A>=t)=1;
    As(:,:,i)=A;
end


% sufficient stats
constants   = get_constants(As,ClassIDs);
phat        = get_ind_edge_params(As,constants);
SigMat      = run_get_fisher_pvals(As,constants);

% parameters
m=mhat;
s=shat;
[coherent wcounter] = coherent_estimator(SigMat,m,s);

fname='cov_BLSA0317_rates';

%% compute correlations

[I J]=ind2sub(size(As(:,:,1)),coherent);

for i=1:constants.s
    ass=As(:,:,i);
    SS(:,i)=ass(coherent);
end

corrSS=corr(SS');



%% plot 
figure(1), clf
fs=12;
imagesc(corrSS), colorbar
axis('square')
title('Correlation Matrix','fontsize',fs)
ylabel('vertex','fontsize',fs)
xlabel('vertex','fontsize',fs)
set(gca,'fontsize',fs)

% print_fig(['../figs/' fname],[6 2]*1.5)

%%
fname='BLSA0317_coclustered_corr';

clc;
k=2; % # of clusters
[row_clust_idx, col_clust_idx,y_index,x_index]=SpectralCoClustering(corrSS,k);

figure(k); clf; colormap(gray);
imagesc(corrSS(y_index,x_index));
% set(gca,'FontSize',14,'FontWeight','bold');
colorbar
title('Correlation Matrix','fontsize',fs)
xlabel('vertex','fontsize',fs)
ylabel('vertex','fontsize',fs)
axis('square')
set(gca,'XTick',[100:100:400],'YTick',[100:100:400])


if savestuff==1
    print_fig(['../../figs/' fname],[2.3 2]*1.5)
end

%% generate synthetic data


mean0=mean(As(:,:,constants.y0),3);
mean1=mean(As(:,:,constants.y1),3);

totalmean=mean(As,3);
totalmean(totalmean==0)=1/(10*constants.s);
totalmean(totalmean==1)=1-1/(10*constants.s);

E0=totalmean;
E0(coherent)=mean0(coherent);
E0(E0==0)=1/(10*constants.s);
E0(E0==1)=1-1/(10*constants.s);

E1=totalmean;
E1(coherent)=mean1(coherent);
E1(E1==0)=1/(10*constants.s);
E1(E1==1)=1-1/(10*constants.s);

s=constants.s;
s0=constants.s0;
s1=constants.s1;
y0=constants.y0;
y1=constants.y1;

As = nan(V,V,s);
As(:,:,y0)=repmat(E0,[1 1 s0]) > rand(V,V,s0);
As(:,:,y1)=repmat(E0,[1 1 s1]) > rand(V,V,s1);

%% synthetic correlation matrix

clear SS
for i=1:constants.s
    ass=As(:,:,i);
    SS(:,i)=ass(coherent);
end
SS0=find(SS==0);
SS1=find(SS==1);
SS(SS0)=rand(size(SS0))*0.01;
SS(SS1)=rand(size(SS1))*0.01+0.99;

corr_syn=corr(SS');

%% cluster synthetic matrix

k=2; % # of clusters
[row_clust_idx, col_clust_idx,y_index,x_index]=SpectralCoClustering(corr_syn,k);

figure(k); clf; colormap(gray);
imagesc(corr_syn(y_index,x_index));
colorbar


%% kstest to compare the two

[h pval_ks] = kstest2(corrSS(:),corr_syn(:))


%%

[coherogram Coherogram wset] = get_coherograms(As,ClassIDs);
figure(3), clf
fs=12;
imagesc(Coherogram(:,1:wcounter))
colorbar
title('data coherogram','fontsize',fs)
ylabel('vertex','fontsize',fs)
xlabel('log threshold','fontsize',fs)
set(gca,'XTick',5:5:25,'XTickLabel',round(10*log10(wset(5:5:25)))/10)
set(gca,'fontsize',fs)

fname='BLSA0317_coherogram';
% print_fig(['../figs/' fname],[4 2]*1.5)

