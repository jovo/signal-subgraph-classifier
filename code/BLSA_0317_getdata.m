clear, clc, %clf

etc.dataset     = 'BLSA_0317/';
etc.dir_data    = '~/Research/data/MRI/BLSA/';
etc.dir_results = [etc.dir_data 'results/'];
etc.dir_figs    = '~/Research/figs/MRI/BLSA/';
subdir          = [etc.dir_data etc.dataset 'countMtx/'];

files   = dir(subdir);      % load csv files storing adjacency matrices
files(1:2)=[];              %% remove . and ..
csvfile = importdata([etc.dir_data  'IACL-blsa-subjectdata-1110.csv']); % load class labels

n       = length(files);    % number of experiments
V       = 70;               % number of vertices
AdjMats = nan(V,V,n);       % pre-allocate memory for adjacency matrices
ClassIDs= nan(1,n);         % pre-allocate memory for class ids


for i=1:n
    A=importdata([subdir files(i).name]);  
    A(isnan(A))=0;
    AdjMats(:,:,i)=A;
    
    for k=2:length(csvfile.textdata)
        if strcmpi(files(i).name(6:7),csvfile.textdata{k,1}(1:2))
            if strcmpi(csvfile.textdata{k,3},'male')
                ClassIDs(i)=1;
            else
                ClassIDs(i)=0;
            end
        end
    end
    
end

save('~/Research/data/MRI/BLSA/BLSA_0317/BLSA_0317_countMtx','AdjMats','ClassIDs')


