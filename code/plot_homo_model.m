% load homo synthetic data and make fig 1
clear, clc, clf

name='homo';

load(['/Users/jovo/Research/data/MRI/BLSA/synthetic/', name '.mat']);

As=sim.As;
ClassIDs=sim.targs;

siz=size(As);
V=siz(1);

%% make true params fig

x.save=1;
x.figdir='/Users/jovo/Research/figs/MRI/BLSA/';
x.fname=[name '_true'];

plot_model(sim.params,x)

%% make est params fig


inds=randperm(siz(3));
n=1000;
inds=inds(1:n);

constants=get_constants(As(:,:,inds),ClassIDs(inds));
phat = get_ind_edge_params(As(:,:,inds),constants);

x.fname=[name '_est'];
plot_model(phat,x)


%% test that subsampling the data is working

err=nan(5,10);
ss=100:100:1000;

for j=1:length(ss)
    
    disp(['s = ', num2str(ss(j))])
    
    for i=1:5
        
        disp(['trial # ', num2str(i)])
        
        inds=randperm(siz(3));
        n=ss(j);
        inds=inds(1:n);
        constants=get_constants(As(:,:,inds),ClassIDs(inds));
        phat = get_ind_edge_params(As(:,:,inds),constants);
        
        err(i,j)=norm([phat.E0 phat.E1]-[sim.params.E0 sim.params.E1],2);
        
    end
end

figure(1), plot(mean(err,1))

%% set algorithms

SigInds=find(M.E0-M.E1);
s=length(SigInds);


alg(1).name='coherent';
alg(1).star_list=1;    % list of # star-vertices
alg(1).edge_list{1}=s;

alg(2).name='incoherent';
alg(2).star_list=V;    % list of # star-vertices
alg(2).edge_list{1}=s;

alg(3).name='naive bayes';
alg(3).star_list=V;    % list of # star-vertices
alg(3).edge_list{1}=choose(V,2);


%% plot signal subgraph estimates

inds=randperm(siz(3));
n=100;
inds=inds(1:n);

constants=get_constants(As(:,:,inds),ClassIDs(inds));
SigMat  = run_get_fisher_pvals(As(:,:,inds),constants);
coh = coherent_estimator(SigMat,alg(1).star_list,alg(1).edge_list{1});
inc = coherent_estimator(SigMat,alg(2).star_list,alg(2).edge_list{1});

phat.tru=zeros(V);
phat.coh=zeros(V);
phat.inc=zeros(V);

phat.tru(SigInds)=1;
phat.coh(coh)=1;
phat.inc(inc)=1;

x.fname='homo_ss';
plot_sig_subgraph(phat,x)

%% get Lhats

ss=20:20:200;

for j=1:length(ss)
    disp(['s = ', num2str(ss(j))])
    for i=1:20
        disp(['trial # ', num2str(i)])
        inds=randperm(siz(3));
        n=ss(j);
        inds=inds(1:n);
        for t=1:3
            [alg(t).Lhat(i,j) alg(t).coherent(i,j)] = coherent_classifier_loop(...
                As(:,:,inds),ClassIDs(inds),alg(t).star_list,alg(t).edge_list,0);
            
%             alg(t).Lhat1(i,j) = coherent_classifier_loop(...
%                 As(:,:,inds),ClassIDs(inds),alg(t).star_list,alg(t).edge_list,1);

            if t~=3
                alg(t).numcorrect(i,j)=numcorrect(alg(t).coherent{i,j},SigInds);
            end     
        end
    end
end

%% plot

figure(2), clf
fs=12;

subplot(121)
hold all
for t=1:3
    errorbar(ss(1:j),nanmean(alg(t).Lhat1), nanstd(alg(t).Lhat1)./ss(1:j))
end
set(gca,'XTick',ss)
set(gca,'Ylim',[0 0.5],'XLim',[0 max(ss)+5])
ylabel('misclassification rate','fontsize',fs)
xlabel('# of training samples','fontsize',fs)
legend(alg(1).name,alg(2).name,alg(3).name,'Location','Best')


subplot(122)
hold all
for t=1:2
    errorbar(ss,1-mean(alg(t).numcorrect/s), std(alg(t).numcorrect/s)./ss)
end
set(gca,'XTick',ss)
set(gca,'Ylim',[0 1],'XLim',[0 max(ss)+5])
ylabel('fraction missed edges','fontsize',fs)
xlabel('# of training samples','fontsize',fs)
legend(alg(1).name,alg(2).name,alg(3).name,'Location','Best')

