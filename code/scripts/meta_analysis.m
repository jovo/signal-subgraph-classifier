%% load real data
clear, clc
datatype='synthetic';

if strcmp(datatype,'synthetic')
    load('../../data/BLSA0317_Count_Lhats');
    fname='BLSA0317_Count_synthetic';
    
    mean0=mean(As(:,:,constants.y0),3);
    mean1=mean(As(:,:,constants.y1),3);
    
    totalmean=mean(As,3);
    E0=totalmean;
    E0(coherent)=mean0(coherent);
    
    E1=totalmean;
    E1(coherent)=mean1(coherent);
    
    s=constants.s;
    s0=constants.s0;
    s1=constants.s1;
    y0=constants.y0;
    y1=constants.y1;
    
    As = nan(V,V,s);
    As(:,:,y0)=repmat(E0,[1 1 s0]) > rand(V,V,s0);
    As(:,:,y1)=repmat(E0,[1 1 s1]) > rand(V,V,s1);
    
elseif strcmp(datatype,'real')
    load('~/Research/data/MRI/BLSA/BLSA_0317/BLSA_0317_countMtx');
    fname='BLSA0317_Count_Lhats';
    t=200;
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
    
end



savestuff=1;
if savestuff==1, save(['../../data/' fname]); end

%% alg stuff
i=0;

i=i+1;
alg(i).name='naive bayes';
alg(i).edge_list=find(tril(ones(V)-diag(ones(V,1)),-1));

i=i+1;
alg(i).name='incoherent';
alg(i).edge_list=1:nchoosek(V,2)/2; %round(logspace(1,log10(nchoosek(V,2)/2),100));

i=i+1;
alg(i).name='coherent';
alg(i).star_list=1:V/2;
for m=1:length(alg(i).star_list)
    alg(i).edge_list{m}=1:V*alg(i).star_list(m)/2;
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

if savestuff==1, save(['../../data/' fname]); end



%% generate other stats
% because matlab saving doesn't work:
% for i=1:nAlgs, L(i).Lhat=Out(i).Lhat; end

% sufficient stats
constants   = get_constants(As,ClassIDs);
phat        = get_ind_edge_params(As,constants);
SigMat      = run_get_fisher_pvals(As,constants);

% parameters
[mhat Imin]=min(I);
shat = J(Imin);
[coherent wcounter] = coherent_estimator(SigMat,mhat,shat);
[coherogram Coherogram wset] = get_coherograms(As,ClassIDs);

if savestuff==1, save(['../../data/' fname]); end

%% do Lhat vs. n in synthetic data

if strcmp(datatype,'synthetic')
    synthetic_rates
end


%% plot rates


% fname='BLSA0317_Count_Lhats';
% load(['../../data/' fname])

clear h
figure(1), clf
gray1=0.5*[1 1 1];
gray2=0.75*[1 1 1];
sa=0.5;
fs=12;
fs2=12;
xmax=1000;

if isempty(strfind(fname,'synthetic'))
    nrows=3;
else
    nrows=2;
end


% incoherent Lhats
h(1)=subplot(nrows,2,1);
semilogx(Lhats{2},'color','k','linewidth',2)
xlabel('log size of signal subgraph','fontsize',fs2)
ylabel([{'misclassification rate'}],'fontsize',fs2)
title(['incoherent estimator'],'fontsize',fs2)
if isempty(strfind(fname,'synthetic'))
    axis([1 1000 0 max(Lhats{2})])
    set(gca,'YTick',[0 .25 .5],'XTick',[1 10 100 1000])
else
    axis([1 1001 0 1])
    set(gca,'YTick',[0 .25 .5 .75 1],'XTick',[1 10 100 1000])
end
grid on

% coherent Lhats
h(2)=subplot(nrows,2,2);
L3hat=Lhats{3};
Lmax=max(L3hat(:));
L3hat(isnan(L3hat))=Lmax;
imagesc(L3hat(:,1:1000)),

if isempty(strfind(fname,'synthetic'))
    %     L3hat(L3hat>0.5)=0.5;
    colorbar('YLim',[min(L3hat(:)) 0.5],'YTick',[min(L3hat(:)) 0.3:0.1:0.5],'YTickLabel',[0.13 0.3:0.1:0.5])
else
    colorbar('YLim',[min(L3hat(:)) max(L3hat(:))],'YTick',[min(L3hat(:)) 0.3:0.2:0.9],'YTickLabel',[0.18 0.3:0.2:0.9])
end
xlabel('size of signal subgraph','fontsize',fs2)
ylabel('# star-vertices','fontsize',fs2)
title('coherent estimator','fontsize',fs)
colormap('gray')


if isempty(strfind(fname,'synthetic'))
    
    % selected coherent Lhats
    h(3)=subplot(323);
    hold all
    
    goodones=find(Lhats{3}<(min(Lhats{3}(:))+1/(n+1)));
    [I J]=ind2sub(size(L3hat),goodones);
    Iunique=unique(I);
    for i=Iunique
        semilogx(Lhats{3}(i,:)','linewidth',2)
    end
    axis([1 1000 0 max(max(Lhats{3}(Iunique,:)))])
    set(gca,'YTick',[0 0.16 .25 .5])
    
    
    set(gca,'XScale','log')
    xlabel('log size of signal subgraph','fontsize',fs2)
    ylabel([{'misclassification rate'}],'fontsize',fs2)
    title(['some coherent estimators'],'fontsize',fs)
    grid on
    
    % zoomed coherent Lhats
    h(4)=subplot(324);
    Lhatzoom=Lhats{3}(min(I):max(I),min(J):max(J));
    imagesc(Lhatzoom), colorbar
    colorbar('YLim',[min(L3hat(:)) 0.5],'YTick',[min(L3hat(:)) 0.3:0.1:0.5],'YTickLabel',[0.13 0.3:0.1:0.5])
    xtick=[400:100:600]; %min(J):100:max(J);
    ytick=(min(I):3:max(I));
    set(gca,'YTick',ytick-min(I),'YTickLabel',ytick)
    set(gca,'XTick',xtick-min(J),'XTickLabel',xtick)
    
    xlabel('size of signal subgraph','fontsize',fs2)
    ylabel('# star-vertices','fontsize',fs2)
    title('zoomed in coherent estimator','fontsize',fs)
    
    % signal subgraph
    h(5)=subplot(325);
    blank=zeros(V);
    blank(coherent)=1;
    imagesc(blank)
    title('coherent signal subgraph estimate','fontsize',fs)
    ylabel('vertex','fontsize',fs2)
    xlabel('vertex','fontsize',fs2)
    
    
    
    % coherogram
    h(6)=subplot(326);
    imagesc(Coherogram(:,1:wcounter))
    colorbar
    title('coherogram','fontsize',fs)
    xlabel('threshold')
    xtick=30:30:wcounter;
    set(gca,'XTick',xtick,'XTickLabel',round(100*wset(xtick))/100)
    
    set(h,'fontsize',fs)
    if savestuff==1, print_fig(['../../figs/' fname '_results'],[5 5]*1.5); end
    
else
    
    for n=1:3
        keeps(n).correct_edges=Out(n).correct_edges;
        keeps(n).Lhat=Out(n).Lhat;
        
        keeps(n).misedge_rate=1-mean(keeps(n).correct_edges,2)/shat;
        keeps(n).misedge_ste=(keeps(n).misedge_rate).*(1-keeps(n).misedge_rate)/sqrt(nTrials);%std(keeps(n).correct_edges,[],2)/(shat*sqrt(nTrials));
        keeps(n).Lhat_avg=mean(keeps(n).Lhat,2);
        keeps(n).Lhat_ste=keeps(n).Lhat_avg.*(1-keeps(n).Lhat_avg)/sqrt(nTrials);%std(keeps(n).Lhat,[],2)/sqrt(nTrials);
    end
    
    gray1=0.5*[1 1 1];
    gray2=0.75*[1 1 1];
    sa=0.5;
    
    %
    h(3)=subplot(nrows,2,3); hold all
    errorbar(ns,keeps(1).misedge_rate,keeps(1).misedge_ste,'k','linewidth',2)
    errorbar(ns,keeps(2).misedge_rate,keeps(2).misedge_ste,'color',gray1,'linewidth',2)
    axis([0 ns(end) 0 1])
    % set(gca,'YTick',[2:2:s])
    xlabel('# training samples','fontsize',fs)
    ylabel('missed-edge rate','fontsize',fs)
    
    h(4)=subplot(nrows,2,4); hold all
    errorbar(ns,keeps(1).Lhat_avg,keeps(1).Lhat_ste,'k','linewidth',2)
    errorbar(ns,keeps(2).Lhat_avg,keeps(2).Lhat_ste,'color',gray1,'linewidth',2)
    errorbar(ns,keeps(3).Lhat_avg,keeps(3).Lhat_ste,'color',gray2,'linewidth',2)
    axis([0 ns(end) 0 0.5])
    xlabel('# training samples','fontsize',fs)
    ylabel('misclassification rate','fontsize',fs)
    % legend(alg(1).name,alg(2).name, alg(3).name,'fontsize',fs)
    legend('coh','inc','nb'); %,'fontsize',fs,'Location','Best')
    set(gca,'YTick',[0.1:0.1:0.5])
    % set(gca,'YScale','log')
    
    if savestuff==1
        print_fig(['../../figs/' fname '_Lhats'],[5 3]*1.5)
    end
    
end



