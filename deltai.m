function d=deltai(alpha,beta,gamma,delta,t)

deltai=delta;
for i=2:t
    deltai(i)=alpha-(beta*gamma)/deltai(i-1);
end
d=deltai;

