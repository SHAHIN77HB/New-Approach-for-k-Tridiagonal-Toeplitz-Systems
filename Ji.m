function J=Ji(n,i)
I=eye(n-2*i);
J=[I(end,:);I(1:end-1,:)];
