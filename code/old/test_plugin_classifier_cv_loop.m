%% synthetic
clear, clc

% # vertices
V=10;

% signal subgraph
Ess = [1 2 4 6 8 9];
thesem = zeros(V);
thesem(5,Ess) = 1;
thesem(Ess,5) = 1;
Ess=find(thesem);

% parameters
p=0.5;
q=0.1;
m=1;
s=numel(Ess);
E0=p*ones(V);
E1=p*ones(V);
E1(Ess)=q;

% sample data
n=100;
As = nan(V,V,n);
As(:,:,1:n/2)=repmat(E0,[1 1 n/2]) > rand(V,V,n/2);
As(:,:,n/2+1:n)=repmat(E1,[1 1 n/2]) > rand(V,V,n/2);
ys=[zeros(1,n/2) ones(1,n/2)];


a=0;

% coherent stuff
a=a+1;
alg(a).name='coherent';
alg(a).star_list=[1 5];
alg(a).edge_list{1}=1:V;
alg(a).edge_list{2}=1:5*V;
if numel(alg(a).edge_list)~=numel(alg(a).star_list), error('edge_list doesnt work with star_list'); end

% incoherent stuff
a=a+1;
alg(a).name='incoherent';
alg(a).edge_list=1:50;

% nb stuff
a=a+1;
alg(a).name='naive bayes';
alg(a).edge_list=find(ones(V)-diag(ones(V,1)));

%%
% loo classify
Out1 = plugin_classifier_cv_loop(As,ys,alg,'loo');


%% in sample stuff

Out2 = plugin_classifier_cv_loop(As,ys,alg,'InSample');

%% hold-out stuff
n=100;
Atst = nan(V,V,n);
Atst(:,:,1:n/2)=repmat(E0,[1 1 n/2]) > rand(V,V,n/2);
Atst(:,:,n/2+1:n)=repmat(E1,[1 1 n/2]) > rand(V,V,n/2);

ytst=[zeros(1,n/2) ones(1,n/2)];

Out3 = plugin_classifier_cv_loop(As,ys,alg,'HoldOut',Atst,ytst);



%% plot
figure(1), clf
for i=1:numel(Out1)
    subplot(1,numel(Out1),i), cla
    hold all
    plot(Out1(i).Lhat','k','linewidth',2)
    plot(Out2(i).Lhat','b')
    plot(Out3(i).Lhat','r')
    legend('loo','in','ho')
end

