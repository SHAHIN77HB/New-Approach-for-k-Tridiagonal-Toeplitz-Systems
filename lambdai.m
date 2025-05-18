function l=lambdai(alpha,beta,gamma,lambda,t)

lambdai=lambda;
for i=2:t
    lambdai(i)=alpha-(beta*gamma)/lambdai(i-1);
end
l=lambdai;

