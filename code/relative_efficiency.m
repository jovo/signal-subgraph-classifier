%% synthetic
clear, clc

% signal subgraph
m=1;
s=8;
Ess = 2:s+1;

% parameters
p=0.1;
q=0.2;

fname=['RE_s' num2str(s) '_p' num2str(round(p*100)) '_q' num2str(round(q*100))];

%% simulate data for coherogram plots
clear h SigMat coherent incoherent coherogram Coherogram wset


ns=[10 50 100];  %[10:10:100]; %%%%50:50:500;
Vs=[25 50]; %[10 25 50]; %10:10:100;
nTrials=100;

inc_cor = nan(length(ns),length(Vs));
coh_cor = inc_cor;

for i=1:length(ns)
    
    disp(['n = ' num2str(ns(i))])
    
    for j=1:length(Vs)
        
        disp(['V = ' num2str(Vs(j))])
        
        % parameters
        V=Vs(j);
        E0=p*ones(V);
        E1=p*ones(V);
        E1(Ess)=q;
        
        for t=1:nTrials
            
            % training data
            ntrn=ns(i);
            Atrn = nan(V,V,ntrn);
            Atrn(:,:,1:2:ntrn)=repmat(E0,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
            Atrn(:,:,2:2:ntrn)=repmat(E1,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
            for nn=1:ntrn, Atrn(:,:,nn)=tril(Atrn(:,:,nn),-1); end
            ytrn=repmat([0 1],1,ntrn/2);
            
            constants   = get_constants(Atrn,ytrn);     % get constants to ease classification code
            SigMat   = run_get_fisher_pvals(Atrn,constants);
            incoherent = incoherent_estimator(SigMat,s);
            [coherent cocount]= coherent_estimator(SigMat,m,s);
            
            inc_cor(i,j,t) = numcorrect(incoherent,Ess);
            coh_cor(i,j,t) = numcorrect(coherent,Ess);
            
        end
    end
    
end

save(['../data/' fname])


%% get stats

% load(['../data/' fname])

inc_avg=mean(inc_cor,3);
coh_avg=mean(coh_cor,3);

inc_std=std(inc_cor,[],3);
coh_std=std(coh_cor,[],3);
%%

figure(max(Vs)), clf, clear h
grays=linspace(0,0.75,length(Vs));
gray1=0.5*[1 1 1];
sa=0.5;
fs=10;

h(1)= subplot(121); hold all
for j=[1:length(ns)]%1:5:length(Vs)
    errorbar(Vs,inc_avg(j,:),inc_std(j,:)./sqrt(nTrials),'-','linewidth',2,'color',grays(j)*[1 1 1])
    errorbar(Vs,coh_avg(j,:),coh_std(j,:)./sqrt(nTrials),'-.','linewidth',2,'color',grays(j)*[1 1 1])
end
title('edges','fontsize',fs)
ylabel('% correct','fontsize',fs)
xlabel('# vertices','fontsize',fs)
lh=legend('incoherent','coherent');
text(.5,.5,'V=10')
text(.5,.5,'V=10')
axis('tight')

set(lh,'Location','Best','fontsize',fs)

h(2)= subplot(122); hold all
for j=1:length(ns)
    errorbar(Vs,inc_avg(j,:),inc_std(j,:)./sqrt(20),'-','linewidth',2,'color',grays(j)*[1 1 1])
    errorbar(Vs,coh_avg(j,:),coh_std(j,:)./sqrt(20),'-.','linewidth',2,'color',grays(j)*[1 1 1])
end
title('classes','fontsize',fs)
ylabel('% correct','fontsize',fs)
xlabel('# vertices','fontsize',fs)

set(h,'fontsize',fs)


% %%
% 
% h(1)=subplot(121); hold all
% errorbar(Vs,inc_avg,inc_std,'k','linewidth',2)
% errorbar(ns,coh_avg,coh_std,'color',gray1,'linewidth',2)
% axis([0 ns(end) 0 s+0.1])
% set(gca,'YTick',[2:2:s])
% xlabel('# training samples','fontsize',fs)
% ylabel('# identified edges','fontsize',fs)
% 
% h(2)=subplot(122); hold all
% errorbar(ns,mean(Out(1).Lhat,2),std(Out(1).Lhat,[],2)/sqrt(nTrials),'k','linewidth',2)
% errorbar(ns,mean(Out(2).Lhat,2),std(Out(2).Lhat,[],2)/sqrt(nTrials),'color',gray1,'linewidth',2)
% errorbar(ns,mean(Out(3).Lhat,2),std(Out(3).Lhat,[],2)/sqrt(nTrials),'color',gray2,'linewidth',2)
% axis([0 ns(end) 0 0.5])
% xlabel('# training samples','fontsize',fs)
% ylabel('misclassification rate','fontsize',fs)
% % legend(alg(1).name,alg(2).name, alg(3).name,'fontsize',fs)
% legend('coh','inc','nb','fontsize',fs,'Location','Best')
% set(gca,'YTick',[0.1:0.1:0.5])
% 
% 
% 
% for i=1:numel(h)
%     op=get(h(i),'OuterPosition');
%     op(4)=op(4)+op(2);
%     op(2)=0;
%     set(h(i),'OuterPosition',op)
%     set(h(i),'fontsize',fs)
% end
% 
% print_fig(['../figs/' fname '_Lhats'],[6 2]*1.5)
% 
% 
