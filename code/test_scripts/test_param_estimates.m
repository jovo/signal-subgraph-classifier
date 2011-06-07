% test different estimators
clear, clc

p=0.01;
V=1000;
n=50;
E0=p*ones(V);
E1=p*ones(V);

As(:,:,1:n/2)=repmat(E0,[1 1 n/2]) > rand(V,V,n/2);
As(:,:,n/2+1:n)=repmat(E1,[1 1 n/2]) > rand(V,V,n/2);
ys=[zeros(1,n/2) ones(1,n/2)];

constants = get_constants(As,ys);

Es=nan(3,2*V^2)';
for i=1:3

    if i==1, type='hack';
    elseif i==2, type='robust';
    elseif i==3, type='map';
    end
    
    temp = get_ind_edge_params(As,constants,type);

    Es(:,i)=[temp.E0(:); temp.E1(:)];

    tit(i)=numel(find(Es(:,i)==0)>0);
    if numel(find(Es(:,i)==0)>0), display(type); end
end
x=[0 linspace(min(Es(:)),max(Es(:)),3)];
hist(Es,x)
legend([num2str(tit(1))],[num2str(tit(2))],[num2str(tit(3))])
ass=sort(Es);
ass(1,:)


%%
t=1000000;
a=1.01;
b=1.01;
alpha=a*ones(t,1);
beta=b*ones(t,1);
R=betarnd(alpha,beta);
hist(R,50)