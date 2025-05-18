function [x,time]=BDk_TriD(T,n,k,b)
R1=[];
R2=zeros(k,floor(n/2));
tic
for i=0:k-1
    r=[];
    for j=1:n
        if rem(j,k)==i
            r=[r j];
        end
    end
    R1=[R1 r];
    R2(i+1,1:length(r))=r;
end

P=k_Tri_permutation(R1);
A=P'*T*P;
f=P'*b;
x=[];
t=0;
for j=1:k
    
    D=A(t+1:t+sum(R2(j,:)>0),t+1:t+sum(R2(j,:)>0));
    c=f(t+1:t+sum(R2(j,:)>0));
    if rem(length(D),2)==0
        if length(D)==2
            [x1,~]=gmres(D,c);
        else
            x1=Evencase(length(D),D(1,1),D(1,2),D(2,1),D(1,1),D(1,1),c);
        end
    else
        if length(D)==1
            x1=c/D;
        else
            x1=Oddcase(length(D),D(1,1),D(1,2),D(2,1),D(1,1),D(1,1),c);
        end
    end
    x=[x;x1];
    t=t+sum(R2(j,:)>0);
end
time=toc;
x=P*x;


