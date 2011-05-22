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


fname=['homo_V' num2str(V) '_s' num2str(s) '_p' num2str(round(p*100)) '_q' num2str(round(q*100))];



%% simulate data for coherogram plots
clear h SigMat coherent incoherent coherogram Coherogram wset

ns=[8 16 64 256]; %[20 50 100 300];
ncols=length(ns);


for i=1:ncols
    
    % training data
    ntrn=ns(i);
    Atrn = nan(V,V,ntrn);
    Atrn(:,:,1:2:ntrn)=repmat(E0,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
    Atrn(:,:,2:2:ntrn)=repmat(E1,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
    for nn=1:ntrn, Atrn(:,:,nn)=tril(Atrn(:,:,nn),-1); end
    
    ytrn=repmat([0 1],1,ntrn/2);

    constants   = get_constants(Atrn,ytrn);     % get constants to ease classification code
    SigMat{i}   = run_get_fisher_pvals(Atrn,constants);
    incoherent{i} = incoherent_estimator(SigMat{i},s);
    [coherent{i} cocount{i}]= coherent_estimator(SigMat{i},m,s);    
    [coherogram{i} Coherogram{i} wset{i}] = get_coherograms(Atrn,ytrn);
    
end

save(['../data/' fname])


%% plot SigMat and coheregrams
% load(['../data/' fname])

figure(4), clf
wh = [.21 .18];
side = 0.07;
bottom=0.07;
rowpos=0.22;
colpos=0.21;
fs=6;
nrows=4;


for i=1:ncols
    
    
    h(1+(i-1)*nrows) = subplot('Position',[0*colpos+side (4-i)*rowpos+bottom wh]);
    
    imagesc(-log10(SigMat{i}))
    axis('square')
    ylabel(['n=', num2str(ns(i))],'fontsize',fs)
    if i<ncols,
        set(gca,'Xtick',[],'YTick',[]),
    else
        xlabel('vertex','fontsize',fs)
        ylabel([{['n=', num2str(ns(i))]}; {'vertex'}],'fontsize',fs)
    end
    if i==1, title([{'significance'};  {'matrix'}],'fontsize',fs), end
    
    h(2+(i-1)*nrows) = subplot('Position',[1*colpos+side (4-i)*rowpos+bottom wh]);
    ass=zeros(V); ass(incoherent{i})=1;
    imagesc(ass)
    axis('square')
    if i<ncols, set(gca,'Xtick',[],'YTick',[]),
    else set(gca,'Ytick',[]), xlabel('vertex','fontsize',fs)
    end
    if i==1, title([{'incoherent'}; {'estimate'}],'fontsize',fs), end
    ylabel(['# correct = ', num2str(numcorrect(incoherent{i},Ess))],'fontsize',fs)
    
    h(3+(i-1)*nrows) = subplot('Position',[2*colpos+side (4-i)*rowpos+bottom wh]);
    ass=zeros(V); ass(coherent{i})=1;
    imagesc(ass)
    axis('square')
    if i<ncols, set(gca,'Xtick',[],'YTick',[]),
    else set(gca,'Ytick',[]), xlabel('vertex','fontsize',fs)
    end
    if i==1, title([{'coherent'}; {'estimate'}],'fontsize',fs), end
    ylabel(['# correct = ', num2str(numcorrect(coherent{i},Ess))],'fontsize',fs)
            
    h(4+(i-1)*nrows) = subplot('Position',[3*colpos+side (4-i)*rowpos+bottom wh]);
    imagesc(Coherogram{i}(:,1:cocount{i}))
    axis('square')
    set(gca,'Ytick',[])
    if i==1, title([{'cumulative'}; {'coherogram'}],'fontsize',fs), end
    if i==ncols, xlabel('log significance','fontsize',fs), end
    if cocount{i}>2,
    set(gca,'XTick',ceil(linspace(2,cocount{i},3)),'XTickLabel',...
        0.1*round(10*log10((wset{i}(ceil(linspace(2,cocount{i},3))))))),
    end
    
end
colormap('gray')
for i=1:length(h), set(h(i),'fontsize',fs), end

print_fig(['../figs/' fname '_SigIncCohErogram'],[5 4])


