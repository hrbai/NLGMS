clear all;
clc;
load('COIL20.mat','X','Y');
%X = double(fea);
%Y = gnd;
[n,d]=size(X);
X=mapminmax(X',0,1);
class_num = length(unique(Y));%�����
dataset = 'COIL20.mat';
k=300;
acc_mean1=0;
nmi_mean1=0;
c=5;%the number of nearest neighbors.
MAXu=10^(10);
epsilon=10^(-3);
data=[];
%% baseline: using all the features
fpath = ['aefs_result_',dataset,'.txt'];
[nmi_mean0,nmi_std0,acc_mean0,acc_std0]=write_baseline_result(X',class_num,Y,fpath);


%mȡֵ(5Ϊ���)
maxF=floor(k);
minF=floor(k/3);
Nk=floor((maxF-minF)/10);
 for nk=Nk
    m=minF+nk*10;
    m = 250;
%lambda1��lambda2��ѭ��(7*7)
for r=2
    lambda1=10^(r)
    for s=2
        lambda2=10^(s)
 tic
 %% ��ֵ

Z=zeros(n,n);
J=zeros(n,n);
E=sparse(d,n);
Y1=zeros(d,n);
Y2=zeros(n,n);
miu=10^(-6);%increase
rgo=1.1;%increase
result=0;
alpha=zeros(1,n);
results=[];
 %S�ĳ�ֵ
distX = L2_distance_1(X,X);
[distX1, idx] = sort(distX,2);
S = zeros(n);
rr = zeros(n,1);
    for i = 1:n
    di = distX1(i,2:c+2);
    rr(i) = 0.5*(c*di(c+1)-sum(di(1:c)));
    id = idx(i,2:c+2);
    S(i,id) = (di(c+1)-di)/(c*di(c+1)-sum(di(1:c))+eps);%(LZR)initialize S 
    end
 W = InitializationW(X,d,m,S,E);
 
for t=1:30
%����J
J=svso(lambda1/miu,Z+Y2/miu);%����ֵ����������J

%Update Z
I=eye(n);
M=X'*X;
Z=(inv(I+M))*(M-X'*E+J+(X'*Y1-Y2)/miu);

% ����W
[W,Ik] = UpdateW(W,X,n,d,k,m,E,S);

%����E
S0 = (S+S')/2;
D0 = diag(sum(S0));
Ls_= D0 - S0;%(LZR)compute Laplacian matrix
I_=eye(n);
NI=4*Ls_+2*lambda2*I_+miu*I_;
E=pinv(W')*(4*W'*X*Ls_+miu*W'*(X-X*Z)+W'*Y1)/NI;

%����S

dij=L2_distance_1(W'*(X-E),W'*(X-E));
[distX1, idx] = sort(distX,2);
S = zeros(n);
rr = zeros(n,1);
    for i = 1:n
    di = distX1(i,2:c+2);
    rr(i) = 0.5*(c*di(c+1)-sum(di(1:c)));
    id = idx(i,2:c+2);
    S(i,id) = (di(c+1)-di)/(c*di(c+1)-sum(di(1:c))+eps);%(LZR)initialize S 
    end
alpha_=sum(rr)/n;
S0 = (S+S')/2;
D0 = diag(sum(S0));
Ls = D0 - S0;%(LZR)compute Laplacian matrix
%����Y1��Y2
Y1=Y1+miu*(X-X*Z-E);
Y2=Y2+miu*(Z-J);
%����miu
miu=min(rgo*miu,MAXu);
%�ж�������������


    
    
    %Ŀ�꺯����ֵ
    result0=result;
    result=2*trace(W'*(X-E)*Ls*(X-E)'*W)+alpha_*sum(sum(S.^2))+lambda1*sum(svd(Z))+lambda2*trace(E'*E)
    Q=abs(result-result0);
    results = [results,Q];
    if Q<epsilon 
            la1=lambda1;
            la2=lambda2;
            out0=[num2str(la1),num2str(la2)];
            disp(out0)
    break
    
    end
    
 end 
 toc
    %% evaluation
XK=X([Ik'],:);
[nmi_mean,nmi_std,acc_mean,acc_std]=write_baseline_result(XK',class_num,Y,fpath);
[nmi_max,acc_max] = maxmin(X',class_num,Y,20);

if acc_mean>acc_mean1
    max=acc_mean;
    acc_mean1=acc_mean;
end

if nmi_mean>nmi_mean1
    max2=nmi_mean;
    nmi_mean1=nmi_mean;
end
data=[data;k;m;lambda1;lambda2;acc_mean;nmi_mean]
    end %lambda2��ѭ��

end %lambda1��ѭ��

 end %k��Ӧ������m��ѭ��
max
max2
%��¼����
fileID=fopen('data1.txt','a+');                 %���ļ� �����ļ�id��

fprintf(fileID,'k m lambda1 lambda2 acc_mean nmi_mean\n');                    %����ַ�

fprintf(fileID,'%d %d %.5f %.5f %.5f %.5f\n',data);             %�����������data

fclose(fileID);  


