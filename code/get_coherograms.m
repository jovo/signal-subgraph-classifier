function [coherogram Coherogram wset] = get_coherograms(As,targs)

constants = get_constants(As,targs);     % get constants to ease classification code
pvals = run_get_fisher_pvals(As,constants);

wset=unique([0; sort(pvals(:))]);
wset(wset>1-1e-3)=[];
coherogram=zeros(constants.n,length(wset));
Coherogram=zeros(constants.n,length(wset));
deg1=zeros(1,constants.n);
for i=1:length(wset)
    w=wset(i);    
    blank=pvals;
    blank(blank>w)=1;
    blank(blank<=w)=0;
    deg0=deg1;
    deg1=constants.n*2-(sum(blank)+sum(blank,2)');
    degdiff=deg1-deg0;
    Coherogram(:,i)=deg1;
    coherogram(:,i)=degdiff;
end