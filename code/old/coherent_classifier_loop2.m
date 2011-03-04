function [coherent incorrects] = coherent_classifier_loop2(star_list,edge_list,SigMat,Atst,phat,ys)


for m=1:length(star_list)
            
            disp(['# star-vertics: ' num2str(star_list(m))])
            
            for s=1:length(edge_list{m})
                
                if mod(s,100)==0, disp(['# edges: ' num2str(edge_list{m}(s))]); end
                
                [coherent{m,s,i} incorrects(m,s,i)] = ...
                    coherent_classifier(SigMat, Atst, phat, star_list(m), edge_list{m}(s), ys(i));
                
            end
        end
 