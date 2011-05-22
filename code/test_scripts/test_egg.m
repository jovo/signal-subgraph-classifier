
%%
clc
ass=reshape(linspace(0.06,1,V^2),[V V]);
ass=rand(V)+0.06;
lt=find(triu(ass,0));
ass(lt)=1;
[foo ix]=sort(ass(:));

thesem=[2 4 5 8];
m=numel(thesem);
s=5;
ass(thesem,thesem)=rand(m,m)*0.06;
ass(lt)=1;
ass(14)=ass(16);

notm=[3 6]; % 7 9 10];
ass(notm)=rand(length(notm),1)*0.06;

imagesc(ass), colorbar

[egghat wcounter] = egg2_estimator(ass,m,s);

[egg sort(egghat)]

length(intersect(egg,egghat))


%%
% clear, clc
% V=10;
% thesem=[2 4 5 8];
% m=numel(thesem);
% s=5;
% 
% ass=rand(V)+0.06;
% ass(ass>1)=1;
% ass(thesem,thesem)=rand(m,m)*0.06;
% lt=find(triu(ass,+1));
% temp=ones(V);
% temp(lt)=ass(lt);
% ass=temp';
% 
% egg=find(ass<=0.06);
% 
% % add noise
% notm=[3 6]; % 7 9 10];
% ass(notm)=rand(length(notm),1)*0.06;
% 
% 
% [egghat wcounter] = egg2_estimator(ass,m,s);

% egg
% sort(egghat)

% imagesc(ass)

