function [NMI,ARI,grps,similarity,Z] = CBLRR(in_X,true_labs,n_space,alpha,beta,mu,a)
label=true_labs;
[X,] = FilterGenesZero(in_X);
X = log10(X +1);
X = mapminmax(X', 0, 1);
NMI=0;
ARI=0;

for i=1:15
     lambda=1+(i-1);
     Z = LRR(X,alpha,lambda,mu,a,beta);
     Z=normFun(Z);
     [localX, ~] = localize(abs(Z));
     min_neighbour = min(sum(localX~=0,1));
     if (min_neighbour>1)          
          similarity = localX+localX';    
          grps = SpectralClustering(similarity, n_space);
          NMI=Cal_NMI(label, grps);
          ARI=Cal_ARI(label, grps);
     end
end
end

function [ result ] = normFun( M )
%normFun: Laplacian normalization
num = size(M,1);
nM = zeros(num,num);
result = zeros(num,num);
    
for i = 1:num
      nM(i,i) = sum(M(i,:));
end

for i = 1:num
       rsum = nM(i,i);
       for j = 1:num
            csum = nM(j,j);
            if((rsum==0)||(csum==0))
                result(i,j) = 0;
            else
                result(i,j) = M(i,j)/sqrt(rsum*csum);
            end
        end
end
end

function [localX,coverage] = localize( C )
%C is the coefficient matrix
%[tmp,ind] = sort(C,1,'descend');
[m,n]=size(C);
localX=C;
coverage=zeros(1,n);
for i=1:n
   thr=C(i,i)/1.5;
   localX(localX(:,i)<thr,i)=0;
   coverage(1,i)=mean(C(i,i)./localX(localX(:,i)>thr,i));
end
end
