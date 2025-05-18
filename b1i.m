function b1 = b1i (b,beta,gamma,deltai,lambdai,i)
if i==1
    b1=[b(end);b(1)];
end
n=length(b);
s=0;
t=0;
 for k=1:i-1
     p=1;
     q=1;
     for j=1:k
         p=p*deltai(i-j);
         q=q*lambdai(i-j);
     end
     s=s+(-1)^k * (beta^k * b(n-(i-1)+k))/p;
     t=t+(-1)^k * (gamma^k *b(i-k))/q;
 end
 b1(1,1)=b(n-(i-1))+s;
 b1(2,1)=b(i)+t;
