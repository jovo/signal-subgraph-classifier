ns0=10;
nTrials=4;
for i=1:length(ns)
    
    ntrn=ns0;
    conv=0;
    while conv==0
        
        for t=1:nTrials
            
            % training data
            Atrn = nan(V,V,ntrn);
            Atrn(:,:,1:2:ntrn)=repmat(E0,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
            Atrn(:,:,2:2:ntrn)=repmat(E1,[1 1 ntrn/2]) > rand(V,V,ntrn/2);
            for nn=1:ntrn, Atrn(:,:,nn)=tril(Atrn(:,:,nn),-1); end
            ytrn=repmat([0 1],1,ntrn/2);
            
            constants   = get_constants(Atrn,ytrn);     % get constants to ease classification code
            SigMat   = run_get_fisher_pvals(Atrn,constants);
            [coherent cocount]= coherent_estimator(SigMat,m,s);
            
            coh(t) = length(intersect(coherent,Ess));
            
        end %nTrials
        
        cohmean=1-mean(coh)/s;
        
        if cohmean<inc_avg(i)
            re(i)=ntrn/ns(i);
            conv=1;
        else
            ntrn=round(ntrn+10);
        end
        
    end %while
    
    
end
