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

savestuff=1;

rootname='_demo_4x4';
fname=[rootname 'homo_V' num2str(V) '_s' num2str(s) '_p' num2str(round(p*100)) '_q' num2str(round(q*100))];


%% simulate data for coherogram plots

ns=[4 16 64 256]; %[20 50 100 300];
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
    [coherent{i} wset{i} wcounter{i}]= coherent_estimator(SigMat{i},m,s);
    [coherogram{i} Coherogram{i}] = get_coherograms(Atrn,ytrn);
    
end


if savestuff, save(['../../data/sim/' fname]), end
    

%% plot SigMat and coheregrams
% load(['../../data/sim/' fname])

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
%         ylabel([{['n=', num2str(ns(i))]}; {'vertex'}],'fontsize',fs)
        ylabel([{['n=', num2str(ns(i))]}],'fontsize',fs)
        set(gca,'YTick',[]),
    end
    if i==1, title([{'negative log'};  {'significance matrix'}],'fontsize',fs), end
    
    h(2+(i-1)*nrows) = subplot('Position',[1*colpos+side (4-i)*rowpos+bottom wh]);
    ass=zeros(V); ass(incoherent{i})=1;
    imagesc(ass)
    axis('square')
    if i<ncols, set(gca,'Xtick',[],'YTick',[]),
    else set(gca,'Ytick',[]), xlabel('vertex','fontsize',fs)
    end
    if i==1, title([{'incoherent'}; {'estimate'}],'fontsize',fs), end
    ylabel(['# correct = ', num2str(length(intersect(incoherent{i},Ess)))],'fontsize',fs)
    
    h(3+(i-1)*nrows) = subplot('Position',[2*colpos+side (4-i)*rowpos+bottom wh]);
    ass=zeros(V); ass(coherent{i})=1;
    imagesc(ass)
    if isempty(coherent{i}), image(ass); end
    axis('square')
    if i<ncols, set(gca,'Xtick',[],'YTick',[]),
    else set(gca,'Ytick',[]), xlabel('vertex','fontsize',fs)
    end
    if i==1, title([{'coherent'}; {'estimate'}],'fontsize',fs), end
    ylabel(['# correct = ', num2str(length(intersect(coherent{i},Ess)))],'fontsize',fs)
    
    h(4+(i-1)*nrows) = subplot('Position',[3*colpos+side (4-i)*rowpos+bottom wh]);
    imagesc(Coherogram{i}(:,1:wcounter{i}))
    axis('square')
    set(gca,'Ytick',[])
    if i==1, title([{'coherogram'}],'fontsize',fs), end
    if i==ncols, xlabel('log significance','fontsize',fs), end
    if wcounter{i}>2,
        set(gca,'XTick',ceil(linspace(2,wcounter{i},3)),'XTickLabel',...
            0.1*round(10*log10((wset{i}(ceil(linspace(2,wcounter{i},3))))))),
    end
    
end
colormap('gray')
for i=1:length(h), set(h(i),'fontsize',fs), end

if savestuff==1
    print_fig(['../../figs/' fname],[5 4])
end
