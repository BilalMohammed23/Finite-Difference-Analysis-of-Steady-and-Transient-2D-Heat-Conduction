function [T,time,gauss_iter]= steady_gauss(n1,n2,k,T_old,T,dx,dy,X,Y,error,tol)
tic                                                                        %start of simulation time
gauss_iter=1;
while(error>tol)                                                           %Convergence loop
    for i=2:(n1-1)
        for j=2:(n2-1)
            T(i,j)=((1/(k* dx^2))*(T(i-1,j)+T(i+1,j)))+ ((1/(k* dy^2))*(T(i,j+1)+T(i,j-1)));
            %T(i,j)=(k1*(T_old(i-1,j)+T_old(i+1,j)))+ (k2*(T_old(i,j+1)+T_old(i,j-1)));   
        end
    end
    error=max(max(abs(T-T_old)));                                          %First we are spilitting the matrix into an array by taking the highest values in each of the column and then again max of this array gives the max error value
    gauss_iter=gauss_iter+1;
    T_old=T;
    figure(1)   
    [c,d]=contourf(X,Y,T);
    colormap('jet')                                                        %plotting the colormap
    colorbar
    title_text = sprintf('Iteration number=%d for Steady state condition using Gauss-Seidel iterative method',gauss_iter);
    title(title_text);
    clabel(c,d)
    xlabel('X')
    ylabel('Y')
    pause(0.03)
end
toc;
time=toc;
end