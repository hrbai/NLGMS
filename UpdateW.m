function [W,Ik] = UpdateW(W,X,n,d,k,m,E,S)



%W = InitializationW(X,d,m,S,E);

  
%求Ls

S0 = (S+S')/2;
D0 = diag(sum(S0));
Ls = D0 - S0;

%求A
[lambda,l]=max(eig(X*Ls*X'));
I=eye(d);
A=lambda*I-(X-E)*Ls*(X-E)';
for t=1:10
%求P
C=pinv(W'*A*W);
P=A*W*C*W'*A;
%求A波浪(AI)
%利用sort函数给数列从小到大排列，找前几个�?��的�?
%[b,i]=sort(a)b为从小到大的数字，i为对应位置�?
p=diag(P);
[b,q]=sort(p,'descend');
% b=wrev(b);
% q=wrev(q);%wrev：�?叙排�?
% q=q(1:k);
% [Ik,z]=sort(q,'descend');
Ik=q(1:k);
AI=A([q(1:k)'],[q(1:k)']);
%求V
[V,D] = eig(AI);
[~,id]=sort(diag(D),'descend');
%  V = fliplr(V);
V=V(:,id(1:m));
%求U
 for row=1:d
    for column=1:k
     if row==q(column)
         U(row,column)=1;
     else U(row,column)=0;
     end
    end
 end

W0=W;
W=U*V;
error=sum(sum((W-W0).^2));
  if error<10^(-3)
     break

   
  end
end
end
