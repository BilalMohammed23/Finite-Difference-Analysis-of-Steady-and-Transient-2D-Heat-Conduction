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
[X,Y] = meshgrid(x,y);

%%
%Values for t and its allocations:
T=ones(n1:n2);                                                             %Forming the initial values with a constant value = 1
T(1,:)=900;                                                                %Bottom Boundary Temperature
T(n1,:)=600;                                                               %Top Boundary Temperature
T(:,1)=400;                                                                %Left Boundary Temperature
T(:,n2)=800;                                                               %Right Boundary Temperature
T(1,n1)=(800+900)/2;                                                       %Average Temperature values at the corresponding nodes
T(1,1)=(400+900)/2;
T(n2,1)=(400+600)/2;
T(n2,n1)=(800+600)/2;
T_old=T;                                                                   %Indicating T values to T_old
T_jac=T;                                                                   %Indicating T values to T_jac which will used in SOR over Jacobian.
T_gauss=T;                                                                 %Indicating T values to T_gauss which will used in SOR over Gauss-Seidel.

%%
%Essential data for the iterative solvers
error=9e9;
tol=1e-4;
alpha=1.2;
k=2*((dx^2)+(dy^2))/(dx^2 * dy^2);                                         %(Value of k for Implicit forms of Steady-state-conduction equation)
%k1=(dy^2)/2*((dx^2)+(dy^2));                                              %(Aliter for above coefficient k)k1=1/(k.dx2)
%k2=(dx^2)/2*((dx^2)+(dy^2));                                              %(Aliter for above coefficient k)k1=1/(k.dy2)

%%
%Solution providing different iterative methods:
iterative=input('Enter the type of solver for steady state conduction problem, 1-Jac, 2-Gauss, 3-SOR_Gauss, 4-SOR_Jac : ')


if iterative== 1                                                           %Jacobian Iterative solver
    [T_steady_jac,time_steady_jac,jacob_iter]= steady_jac(n1,n2,k,T_old,T,dx,dy,X,Y,error,tol)
    fprintf('Time taken by the Jacobian Iterative solver to solve = %d seconds.',time_steady_jac);
    fprintf('No. of iterations done to solve the 2D Conduction equation = %d',jacob_iter);
end


if iterative== 2                                                           %Gauss-Seidel Iterative solver
    [T_steady_gauss,time_steady_gauss,gauss_iter]= steady_gauss(n1,n2,k,T_old,T,dx,dy,X,Y,error,tol)
    fprintf('Time taken by the Gauss-Seidel Iterative solver to solve = %d seconds.',time_steady_gauss);
    fprintf('No. of iterations done to solve the 2D Conduction equation = %d',gauss_iter);
end


if iterative== 3                                                           %Successive-Over-Relaxation over GS iterative solver
    [T_steady_sor_gauss,time_steady_sor_gauss,sor_gauss_iter]= steady_sor_gauss(n1,n2,k,T_old,T,dx,dy,X,Y,error,tol,alpha,T_gauss)
    fprintf('Time taken by the SOR_Gauss-Seidel Iterative solver to solve = %d seconds.',time_steady_sor_gauss);
    fprintf('No. of iterations done to solve the 2D Conduction equation = %d',sor_gauss_iter);
end


%Unstable for any values above 1 for alpha
if iterative==4                                                            %Successive-Over-Relaxation over Jacobian iterative solver
    [T_steady_sor_jac,time_steady_sor_jac,sor_jac_iter]= steady_sor_jac(n1,n2,k,T_old,T,dx,dy,X,Y,error,tol,alpha,T_jac)
    fprintf('Time taken by the SOR_Jacobian Iterative solver to solve = %d seconds.',time_steady_sor_jac);
    fprintf('No. of iterations done to solve the 2D Conduction equation = %d',sor_jac_iter);
end
