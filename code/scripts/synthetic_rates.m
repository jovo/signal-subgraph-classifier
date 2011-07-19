clear, clc
load('../../data/results/BLSA0317_Count_Lhats');

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

fname='synthetic_BLSA0317_rates';

%% generate data, classify, and get stats
clear alg out Outtemp ns h
% alg stuff
alg(1).name='coherent';
alg(1).star_list=mhat;
alg(1).edge_list{1}=shat;

alg(2).name='incoherent';
alg(2).edge_list=shat;

alg(3).name='naive bayes';
alg(3).edge_list=find(ones(V)-diag(ones(V,1)));

nAlgs=numel(alg);
nTrials=20;
ns=[10:10:100];


for i=1:length(ns)
    
    disp(['ntrn ', num2str(ns(i))])
    
    for j=1:nTrials
        
        disp(['trial ', num2str(j)])
        % training data
        ntrn=ns(i);
        Atrn = nan(V,V,ntrn);
        Atrn(:,:,1:2:ntrn)=repmat(E0,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
        Atrn(:,:,2:2:ntrn)=repmat(E1,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
        
        ytrn=repmat([0 1],1,ntrn/2);
        
        % testing data
        ntst=100;
        Atst = nan(V,V,ntst);
        Atst(:,:,1:2:ntst)=repmat(E0,[1 1 ntst/2]) > rand(V,V,ntst/2);
        Atst(:,:,2:2:ntst)=repmat(E1,[1 1 ntst/2]) > rand(V,V,ntst/2);
        ytst=repmat([0 1],1,ntst/2);
        
        % classify
        Outtemp = plugin_classifier_cv_loop(Atrn,ytrn,alg,'HoldOut',Atst,ytst);
        
        % get stats
        Out(1).Lhat(i,j)   = Outtemp(1).Lhat;
        Out(1).correct_edges(i,j) = length(intersect(Outtemp(1).subspace{1},coherent));
        
        Out(2).Lhat(i,j)   = Outtemp(2).Lhat;
        Out(2).correct_edges(i,j) = length(intersect(Outtemp(2).subspace{1},coherent));
        
        Out(3).Lhat(i,j)   = Outtemp(3).Lhat;
    end
    
end

% for i=1:nAlgs
%     Lhats{i}=Out(i).Lhat;
% end

if savestuff==1, save(['../../data/' fname]); end


%% plot rates

% load(['../data/' fname])

% figure(1), clf, hold all
% 
% for n=1:3
%     keeps(n).correct_edges=Out(n).correct_edges;
%     keeps(n).Lhat=Out(n).Lhat;
%     
%     keeps(n).misedge_rate=1-mean(keeps(n).correct_edges,2)/shat;
%     keeps(n).misedge_ste=(keeps(n).misedge_rate).*(1-keeps(n).misedge_rate)/sqrt(nTrials);%std(keeps(n).correct_edges,[],2)/(shat*sqrt(nTrials));
%     keeps(n).Lhat_avg=mean(keeps(n).Lhat,2);
%     keeps(n).Lhat_ste=keeps(n).Lhat_avg.*(1-keeps(n).Lhat_avg)/sqrt(nTrials);%std(keeps(n).Lhat,[],2)/sqrt(nTrials);
% end
% 
% 
% 
% gray1=0.5*[1 1 1];
% gray2=0.75*[1 1 1];
% sa=0.5;
% fs=12;
% 
% 
% %
% h(1)=subplot(211); hold all
% errorbar(ns,keeps(1).misedge_rate,keeps(1).misedge_ste,'k','linewidth',2)
% errorbar(ns,keeps(2).misedge_rate,keeps(2).misedge_ste,'color',gray1,'linewidth',2)
% axis([0 ns(end) 0 1])
% % set(gca,'YTick',[2:2:s])
% xlabel('# training samples','fontsize',fs)
% ylabel('missed-edge rate','fontsize',fs)
% 
% h(2)=subplot(212); hold all
% errorbar(ns,keeps(1).Lhat_avg,keeps(1).Lhat_ste,'k','linewidth',2)
% errorbar(ns,keeps(2).Lhat_avg,keeps(2).Lhat_ste,'color',gray1,'linewidth',2)
% errorbar(ns,keeps(3).Lhat_avg,keeps(3).Lhat_ste,'color',gray2,'linewidth',2)
% axis([0 ns(end) 0 0.5])
% xlabel('# training samples','fontsize',fs)
% ylabel('misclassification rate','fontsize',fs)
% % legend(alg(1).name,alg(2).name, alg(3).name,'fontsize',fs)
% legend('coh','inc','nb'); %,'fontsize',fs,'Location','Best')
% set(gca,'YTick',[0.1:0.1:0.5])
% % set(gca,'YScale','log')
% 
% if savestuff==1
%     print_fig(['../../figs/' fname '_Lhats'],[3 3]*1.5)
% end
% 
