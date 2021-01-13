function [W,H,err]=ATV_NMF(X,options,k,S,D)  
maxiter =options.maxIter;
[i,u]=size(X);
r = k;                                 
W=rand(i,r);                            
H=rand(r,u);                                 
alpha =options.alpha;
for iter=1:maxiter
   divHH=ATV_div(H);%%        
   W=W.*((X*H')./(W*H*H'));  
    W = max(W,eps);
    H=H.*((W'*X+alpha*H*S+divHH)./(W'*W*H+alpha*H*D));
    H = max(H,eps);   
    err(iter)=norm(X-W*H,'fro')^2;
end
