clear all
close all
clc
%%
%indication of n,x,y terms
n1=input('Enter the number of nodal points(suggestion:enter less than 30 to not exceed the simulation time) : '); %higher the n value, higher is the smoothness and accuracy of the contour;
n2=n1;                                                                     %Indicating the same value for y-axis nodes
x=linspace(0,1,n1);
y=linspace(0,1,n2);
dx=x(2)-x(1);
dy=y(2)-y(1);

%%
%Values for t and its allocations:
T=300*ones(n1:n2);                                                         %Forming the initial values with a constant value = 1
T(1,:)=900;                                                                %Bottom Boundary Temperature
T(n1,:)=600;                                                               %Top Boundary Temperature
T(:,1)=400;                                                                %Left Boundary Temperature
T(:,n2)=800;                                                               %Right Boundary Temperature
T(1,n1)=(800+900)/2;                                                       %Average Temperature values at the corresponding nodes
T(1,1)=(400+900)/2;
T(n2,1)=(400+600)/2;
T(n2,n1)=(800+600)/2;
T_old=T;                                                                   %Indicating T values to T_old
T_initial=T_old;                                                           %Indicating T_old or T values to T_initial which is used for each iteration for time loop.
T_gauss=T;                                                                 %Indicating T values to T_jac which will used in SOR over Gauss-Seidel.
T_jac=T;                                                                   %Indicating T values to T_jac which will used in SOR over Jacobian.

%%
%Essential data for the iterative solvers
error=9e9;
tol=1e-4;                                                                  %absolute error criterion
diff=1.4;                                                                  %diffusion term  %As this value increases we may need increased simulation time for convergence
dt=1e-3;                                                                   %time interval dt
alpha=1.2;                                                                 %Scaling factor
k1= diff * (dt/ dx^2);
k2= diff * (dt/ dy^2);
term1=(1+ 2*k1 + 2*k2)^-1;
term2=k1*term1;
term3=k2*term1;

%%
%Solution providing different iterative methods:
%Implicit forms:
iterative=input('Enter the type of solver for steady state conduction problem, 1-Jac, 2-Gauss, 3-SOR_Gauss, 4-SOR_Jac, 5-Explicit : ')

if iterative==1                                                            %Jacobian Iterative solver
    jac_iter=1
    [T_transient_jac,time_transient_jac,jacob_iter]=transient_jac(n1,n2,T_old,T_initial,T,dx,dy,x,y,tol,term1,term2,term3,jac_iter)
    fprintf('Time taken by the Jacobian Iterative solver to solve = %d seconds.',time_transient_jac);
    fprintf(' No. of iterations done to solve the 2D Conduction equation = %d',jacob_iter);
end

if iterative==2                                                            %Gauss-Seidel Iterative solver
    gauss_iter=1;
    [T_transient_gauss,time_transient_gauss,gauss_iter]=transient_gauss(n1,n2,T_old,T_initial,T,dx,dy,x,y,tol,term1,term2,term3,gauss_iter)
    fprintf('Time taken by the Gauss-Seidel Iterative solver to solve = %d seconds.',time_transient_gauss);
    fprintf(' No. of iterations done to solve the 2D Conduction equation = %d',gauss_iter);
end

if iterative==3                                                            %Successive-Over-Relaxation over GS iterative solver
    sor_gauss_iter=1;
    [T_transient_sor_gauss,time_transient_sor_gauss,sor_gauss_iter]=transient_sor_gauss(n1,n2,T_old,T_initial,T,x,y,tol,term1,term2,term3,alpha,sor_gauss_iter)
    fprintf('Time taken by the SOR_Gauss-Seidel Iterative solver to solve = %d seconds.',time_transient_sor_gauss);
    fprintf(' No. of iterations done to solve the 2D Conduction equation = %d',sor_gauss_iter);
end

if iterative==4                                                            %Successive-Over-Relaxation over Jacobian iterative solver
    sor_jac_iter=1;
    [T_transient_sor_jac,time_transient_sor_jac,sor_jac_iter]=transient_sor_jac(n1,n2,T_old,T_initial,T,x,y,tol,term1,term2,term3,alpha,sor_jac_iter)
    fprintf('Time taken by the SOR_Jacobian Iterative solver to solve = %d seconds.',time_transient_sor_jac);
    fprintf(' No. of iterations done to solve the 2D Conduction equation = %d',sor_jac_iter);
end


if iterative==5                                                            %Explicit form of unsteady state conduction
    exp_iter=1;
    [T_transient_explicit,time_transient_explicit,explicit_iter]=transient_explicit(n1,n2,T_old,T_initial,T,x,y,tol,k1,k2,exp_iter)
    fprintf('Time taken by the Iterative Explicit solver to solve = %d seconds.',time_transient_explicit);
    fprintf(' No. of iterations done to solve the 2D Conduction equation = %d',explicit_iter);
end