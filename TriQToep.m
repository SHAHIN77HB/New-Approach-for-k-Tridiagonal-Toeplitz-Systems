function T=TriQToep(n,alpha,beta,gamma,lambda,delta)
if n==3
    T=[alpha beta 0 ; gamma alpha beta ; 0 gamma alpha];
elseif n==2
    T=[alpha beta ; gamma alpha];
else
    A=toeplitz([alpha gamma zeros(1,n-4)],[alpha beta zeros(1,n-4)]);
    T1=[lambda [beta zeros(1,n-3)] 0];
    T2=[[gamma;zeros(n-3,1)] A [zeros(n-3,1); beta]];
    T3=[0 [zeros(1,n-3) gamma] delta];
    T=[T1;T2;T3];
end