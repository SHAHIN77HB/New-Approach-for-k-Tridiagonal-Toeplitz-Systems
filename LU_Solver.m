function [x,err_lu,time]=LU_Solver(A,b)
tic
[~,l,u]=plu(A);
y=forward_Substitution_System_Solver(l,b);
x=Backward_Substitution_System_Solver_lu(u,y);
time=toc;
err_lu=norm(A*x-b)/norm(b);