function [grps,similarity,Z] = CBLRR(alpha,mu,a,beta,data,label)
n_space = length(unique(label));
    for i=1:15
        lambda=1+(i-1);
        Z = LRR(data,alpha,lambda,mu,a,beta);
        Z=normFun(Z);
        [localX, ~] = localize(abs(Z));
        min_neighbour = min(sum(localX~=0,1));
        if (min_neighbour>1)          
            similarity = localX+localX';    
            grps = SpectralClustering(similarity, n_space);
        end
    end
end