function n = numcorrect(x,y)

n=0;
for i=1:length(x)
    for j=1:length(y)
        if x(i)==y(j), n=n+1; end
    end
end
