%% load real data
clear, clc
datatype='count';

if strcmp(datatype,'FA')
    load('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/BLSA_0317_FAMtx');
    fname='BLSA0317_FA';
    t=0.4;
else
    load('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/BLSA_0317_CountMtx');
    fname='BLSA0317_Count';
    t=200;
end


siz=size(AdjMats);
n=siz(3);           % # experiments
V=siz(1);           % # vertices


As=0*AdjMats;       % threshold
for i=1:n
    A=AdjMats(:,:,i);
    A(A<t)=0;
    A(A>=t)=1;
    As(:,:,i)=tril(A,-1);
end



%% alg stuff
i=0;

i=i+1;
alg(i).name='naive bayes';
alg(i).edge_list=find(tril(ones(V)-diag(ones(V,1)),-1));

i=i+1;
alg(i).name='incoherent';
alg(i).edge_list=1:choose(V,2)/2; %round(logspace(1,log10(choose(V,2)/2),100));

i=i+1;
alg(i).name='coherent';
alg(i).star_list=1:V/2;
for m=1:length(alg(i).star_list)
    alg(i).edge_list{m}=1:V*alg(i).star_list(m);
end

% i=i+1;
% alg(i).name='egg';
% alg(i).star_list=1:V/2;
% for m=1:length(alg(i).star_list)
%     alg(i).edge_list{m}=1:V*alg(i).star_list(m)/2; unique(round(logspace(0,log10(V*alg(i).star_list(m)/3),100))); %
% end

nAlgs=numel(alg);

%% classify
Out = plugin_classifier_cv_loop(As,ClassIDs,alg,'loo');

for i=1:nAlgs
    Lhats{i}=Out(i).Lhat;
end

save(['../../data/' fname])


%% plot rates


% fname='BLSA0317_Count';
% load(['../../data/' fname])

clear h
figure(1), clf

gray1=0.5*[1 1 1];
gray2=0.75*[1 1 1];
sa=0.5;
fs=10;
xmax=1000;
h(1)=subplot(321);
semilogx(Lhats{2},'color','k','linewidth',2)
xlabel('log size of signal subgraph','fontsize',fs)
ylabel([{'misclassification'};  {'rate'}],'fontsize',fs)
title(['incoherent estimator'])
axis([1 1000 0 0.5])
set(gca,'YTick',[0 .25 .5])
grid on

h(2)=subplot(322);
L3hat=Lhats{3};
L3hat(L3hat>0.5)=0.5;
Lmax=max(L3hat(:));
L3hat(isnan(L3hat))=Lmax;
imagesc(L3hat), colorbar
xlabel('size of signal subgraph','fontsize',fs)
ylabel('# star-vertices','fontsize',fs)
title('coherent estimator')


h(3)=subplot(323);
goodones=find(L3hat<(min(L3hat(:))+1/(n+1)));
[I J]=ind2sub(size(L3hat),goodones);
Iunique=unique(I);
hold all
for i=Iunique
    semilogx(L3hat(i,:)','linewidth',2)
end
set(gca,'XScale','log')
xlabel('log size of signal subgraph','fontsize',fs)
ylabel([{'misclassification'};  {'rate'}],'fontsize',fs)
title(['some coherent estimators'])
axis([1 1000 0 0.5])
set(gca,'YTick',[0 .25 .5])
grid on

h(4)=subplot(324);
imagesc(Lhats{3}(min(I):max(I),min(J):max(J))), colorbar
xlabel('size of signal subgraph','fontsize',fs)
ylabel('# star-vertices','fontsize',fs)
title('zoomed in coherent estimator')
colormap('default')
ytick=(min(I):3:max(I));
set(gca,'YTick',ytick-min(I),'YTickLabel',ytick)
xtick=min(J):50:max(J);
set(gca,'XTick',xtick-min(J),'XTickLabel',xtick)


h(5)=subplot(325);
blank=zeros(V);
blank(coherent)=1;
imagesc(blank)
title('coherent signal subgraph estimate','fontsize',fs)

h(6)=subplot(326);
imagesc(Coherogram(:,1:wcounter))
colorbar
title('coherogram','fontsize',fs)
xlabel('threshold')
xtick=30:30:wcounter;
set(gca,'XTick',xtick,'XTickLabel',round(100*wset(xtick))/100)

%
set(h,'fontsize',fs)
print_fig(['../figs/' fname '_results'],[5 3]*1.5)


%% generate other stats
% because matlab saving doesn't work:
for i=1:nAlgs, L(i).Lhat=Out(i).Lhat; end

% sufficient stats
constants   = get_constants(As,ClassIDs);
phat        = get_ind_edge_params(As,constants);
SigMat      = run_get_fisher_pvals(As,constants);

% parameters 
[mhat Imin]=min(I);
shat = J(Imin);
[coherent wcounter] = coherent_estimator(SigMat,mhat,shat);
[coherogram Coherogram wset] = get_coherograms(As,ClassIDs);

save(['../data/' fname])
