clear
clc
%%
%syms x y z
n=2^10;
k=2^10-10;
h=floor(n/k);
alpha=5;beta=1;gamma=2;
T=toeplitz([alpha zeros(1,k-1) gamma zeros(1,n-k-1)],[alpha ;zeros(k-1,1); beta ;zeros(n-k-1,1)]);


%xx=rand(n,1);
%b=T*xx;
b=rand(n,1);
%b=T*kron(ones(n/2,1),[-1;1]);
%% BDk_TriD Tests 

for i=1:10
    [x,time(i)]=BDk_TriD(T,n,k,b);
    err(i)=norm((b-T*x))/norm(b);
    [x_lu,err_lu(i),time_lu(i)]=LU_Solver(T,b);
end
ave_time=sum(time)/10
ave_err=sum(err)/10
ave_time_lu=sum(time_lu)/10
ave_err_lu=sum(err_lu)/10
