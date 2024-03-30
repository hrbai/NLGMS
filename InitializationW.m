function bestW = InitializationW(X,d,m,S,E)




bestofv=0;
bestW=zeros(d,m);

for j=1:10%随机生成初始值，重复10次
  
W=RandOrthMat(d, m, 1);%随机生成满足W'W=I的矩阵W0作为W的初值

%求Ls
S0 = (S+S')/2;
D0 = diag(sum(S0));
Ls = D0 - S0;

%求A
[lambda,l]=max(eig(X*Ls*X'));
I=eye(d);
A=lambda*I-(X-E)*Ls*(X-E)';

ofv=trace(W'*A*W);
    if(ofv>bestofv)
          	bestofv = ofv;
          
            bestW   =  W;
          
    end
end

end










