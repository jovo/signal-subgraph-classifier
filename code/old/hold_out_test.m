% load homo synthetic data and make fig 1
clear, clc, clf

name='homo1';

load(['/Users/jovo/Research/data/MRI/BLSA/synthetic/', name '.mat']);

Atrn=trn.As;
ytrn=trn.targs;

Atst=tst.As;
ytst=tst.targs;

siz=size(Atrn);
V=siz(1);

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


%% get Lhats

ss=500;

for j=1:length(ss)
    disp(['s = ', num2str(ss(j))])
    for i=1:length(ytst)
        disp(['trial # ', num2str(i)])
        inds=randperm(siz(3));
        n=ss(j);
        inds=inds(1:n);
        for t=1:3
            [alg(t).Lhat(i,j) alg(t).coherent{i,j}] = ...
                coherent_classifier_loop(Atrn(:,:,inds),ytrn(inds),alg(t).star_list,alg(t).edge_list,2,Atst,ytst);
            
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

