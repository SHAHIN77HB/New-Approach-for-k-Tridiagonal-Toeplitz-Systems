function [x,err,time]=Evencase(n,alpha,beta,gamma,lambda,delta,b)
t=(n-4)/2+1;
T=TriQToep(n,alpha,beta,gamma,lambda,delta);
J=Ji(n,0);
%%
tic
l=lambda;
d=delta;
l=lambdai(alpha,beta,gamma,lambda,t);
d=deltai(alpha,beta,gamma,delta,t);
At=[alpha beta ; gamma alpha];
Wt=[0 gamma ; beta 0];
Pt=Pi(beta,gamma,n,t);
M=At-Wt*invDi(l,d,t)*Pt;
c=b2i(b,t)-Wt*invDi(l,d,t)*b1i(b,beta,gamma,d,l,t);
[x2,~]=gmres(M,c);
x1=invDi(l,d,t)*( b1i(b,beta,gamma,d,l,t)- Pt*x2);
xt=[x1;x2];
%%
for i=t-1:-1:1
    x2(1:length(xt),t-(i-1))=Ji(n,i)'*xt;
    x1(:,t-(i-1))=invDi(l,d,i)*( b1i(b,beta,gamma,d,l,i)- Pi(beta,gamma,n,i)*x2(:,t-(i-1)));
    xt=[x1(:,end);x2(:,end)];
end
time=toc;
x=J'*xt(:,end);
err=norm(T*x-b)/norm(b);





