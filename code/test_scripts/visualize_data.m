clear, clc


%% load data
% data generated by /Users/jovo/Research/code/necog/analysis/MRI/BLSA/BLSA_0317_getdata
% raw data in /Users/jovo/Research/data/MRI/BLSA/BLSA_0317

load('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/BLSA_0317_FAMtx');
load('/Users/jovo/Research/data/MRI/BLSA/BLSA_0317/BLSA_0317_countMtx');

%% visualize each graph

n=numel(ClassIDs);

for i=1:n
    
    A=AdjMats(:,:,i);
    imagesc(A), colorbar
    
    A(A==0)=[];
    small=find(A<0.2);
    
    large=find(A>1);
    
    title(['i= ', num2str(i), ', <0.2 ', num2str(length(small)), ', >1 ', num2str(length(large))]) 
    keyboard
    
end

%%

As=AdjMats(:);

subplot(121), hist(As,100)

As(As==0)=[];
subplot(122), hist(log(As),100)

% seems bimodal-ish with a local min at around log(As)==2.6, 
% there min at As=exp(2.6) \approx 13.4

%% make plots again using counts but thresholding

clf

for i=1:n
    
    A=AdjMats(:,:,i);

    A(A<exp(2.6))=0;
    A(A>=exp(2.6))=1;
    
    imagesc(A), colorbar
    
    num_edges=sum(A(:)/2);
    
    
    title(['i= ', num2str(i), ', # edges ', num2str(num_edges), ', sparsity= ' num2str(num_edges/2415)]) 
    keyboard
    
end

% seems like about 60% of edges remain when thresholding at 13.4.
% notably, 60% is fairly constant across experiments