%% synthetic
clear, clc


% # neurons
V=70;

% signal subgraph
m=1;
s=20;

Ess = randperm(70);
Ess = Ess(1:s);
thesem = zeros(V);
vstar=round(V/2);
thesem(vstar,Ess) = 1;
thesem(Ess,vstar) = 1;
thesem=tril(thesem,-1);
Ess=find(thesem);

% parameters
p=0.1;
q=0.3;
E0=p*ones(V);
E1=p*ones(V);
E1(Ess)=q;

savestuff=0;

fname=['homo_V' num2str(V) '_s' num2str(s) '_p' num2str(round(p*100)) '_q' num2str(round(q*100))];


%% generate data, classify, and get stats

constraints{1}=NaN;         % use all edges
constraints{2}=19;          % use best 5 edges
constraints{3}=[1 19];     % use best 5 edges incident to 1 vertex


nTrials=100;
ns=[10:10:50 100:50:300];


for i=1:length(ns)
    
    disp(['ntrn ', num2str(ns(i))])
    
    for j=1:nTrials
        
        disp(['trial ', num2str(j)])
        % training data
        ntrn=ns(i);
        Atrn = nan(V,V,ntrn);
        Atrn(:,:,1:2:ntrn)=repmat(E0,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
        Atrn(:,:,2:2:ntrn)=repmat(E1,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
        for nn=1:ntrn, Atrn(:,:,nn)=tril(Atrn(:,:,nn),-1); end
        
        ytrn=repmat([0 1],1,ntrn/2);
        
        % testing data
        ntst=100;
        Atst = nan(V,V,ntst);
        Atst(:,:,1:2:ntst)=repmat(E0,[1 1 ntst/2]) > rand(V,V,ntst/2);
        Atst(:,:,2:2:ntst)=repmat(E1,[1 1 ntst/2]) > rand(V,V,ntst/2);
        for nn=1:ntst, Atst(:,:,nn)=tril(Atst(:,:,nn),-1); end
        ytst=repmat([0 1],1,ntst/2);
        
        
        % signal-subgraph classify
        [Lhat incorrects subspace tElapsed] = xval_SigSub_classifier(Atrn,ytrn,constraints,'HoldOut',Atst,ytst);
        
        
        % logistic regression with L1 penalty classify
        tic
        [B,FitInfo] = lassoglm([ones(1, ntrn);  reshape(Atrn,[V^2 ntrn])]',ytrn','binomial','DFmax',length(Ess));
        yhat=nan(ntst,1);
        parfor k=1:ntst
            yhat(k) = glmval(B(:,1),reshape(Atst(:,:,k),[V^2 1])','logit');
        end
        lassoTime=toc;
        
        % get stats
        Out(1).Lhat(i,j)    = Lhat(1);
        Out(1).time(i,j)    = tElapsed(1);
        
        Out(2).Lhat(i,j)    = Lhat(2);
        Out(2).correct_edges(i,j) = length(intersect(subspace{2},Ess));
        Out(2).time(i,j)    = tElapsed(2);

        Out(3).Lhat(i,j)    = Lhat(3);
        Out(3).correct_edges(i,j) = length(intersect(subspace{3},Ess));
        Out(3).time(i,j)    = tElapsed(3);

        Out(4).Lhat(i,j)    = sum(yhat>0.75~=ytst')/ntst;
        Out(4).correct_edges(i,j) = length(intersect(find(B(:,1))-1,Ess));
        Out(4).time(i,j)    = lassoTime;

    end
    
end

%%

for n=1:4
    keeps(n).correct_edges=Out(n).correct_edges;
    keeps(n).Lhat=Out(n).Lhat;
    
    keeps(n).misedge_rate=1-mean(keeps(n).correct_edges,2)/s;
    keeps(n).misedge_ste=std(keeps(n).correct_edges,[],2)/(s*sqrt(nTrials));

    keeps(n).Lhat_avg=nanmean(keeps(n).Lhat,2);
    keeps(n).Lhat_ste=nanstd(keeps(n).Lhat,[],2)/sqrt(nTrials);

    keeps(n).time_med=nanmedian(Out(n).time,2);
    keeps(n).time_mad=mad(Out(n).time')';
    keeps(n).time_iqr=iqr(Out(n).time')';

end

if savestuff==1
    save(['../../data/' fname])
end

%% plot rates
% load(['../../data/' fname])

figure(1), clf

gray1=0.5*[1 1 1];
gray2=0.75*[1 1 1];
sa=0.5;
fs=12;


%
h(1)=subplot(311); hold all
errorbar(ns,keeps(2).misedge_rate,keeps(2).misedge_ste,'color',gray1,'linewidth',2)
errorbar(ns,keeps(3).misedge_rate,keeps(3).misedge_ste,'color','k','linewidth',2)
errorbar(ns,keeps(4).misedge_rate,keeps(4).misedge_ste,'--','color',gray2,'linewidth',2)
axis([0 ns(end)+5 0 1])
% set(gca,'YTick',[2:2:s])
xlabel('# training samples','fontsize',fs)
ylabel('missed-edge rate','fontsize',fs)

h(2)=subplot(312); hold all
errorbar(ns,keeps(1).Lhat_avg,keeps(1).Lhat_ste,'--','color','k','linewidth',2)
errorbar(ns,keeps(2).Lhat_avg,keeps(2).Lhat_ste,'color',gray1,'linewidth',2)
errorbar(ns,keeps(3).Lhat_avg,keeps(3).Lhat_ste,'color','k','linewidth',2)
errorbar(ns,keeps(4).Lhat_avg,keeps(4).Lhat_ste,'--','color',gray2,'linewidth',2)
axis([0 ns(end)+5 0 0.5])
xlabel('# training samples','fontsize',fs)
ylabel('misclassification rate','fontsize',fs)
set(gca,'YTick',[0.1:0.1:0.5])



h(3)=subplot(313); hold all
errorbar(ns(5:end),keeps(1).time_med(5:end),keeps(1).time_iqr(5:end)/2,'-.','color','k','linewidth',2)
errorbar(ns(5:end),keeps(2).time_med(5:end),keeps(2).time_iqr(5:end)/2,'color',gray1,'linewidth',2)
errorbar(ns(5:end),keeps(3).time_med(5:end),keeps(3).time_iqr(5:end)/2,'color','k','linewidth',2)
errorbar(ns(5:end),keeps(4).time_med(5:end),keeps(4).time_iqr(5:end)/2,'--','color',gray2,'linewidth',2)
axis([ns(1) ns(end)+5 0 max(keeps(4).time_med(5:end)+keeps(4).time_iqr(5:end)/2)])
xlabel('# training samples','fontsize',fs)
ylabel('time (sec)','fontsize',fs)

% legend(alg(1).name,alg(2).name, alg(3).name,'fontsize',fs)
legend('nb', 'inc', 'coh','lasso','Location','NorthWest')
set(h(3),'yscale','log','ytick',[0.1 0.5 1], 'yticklabel',[0.1 0.5 1])

if savestuff==1
    print_fig(['../../figs/fig3_' fname '_Lhats'],[3 3]*2)
end

