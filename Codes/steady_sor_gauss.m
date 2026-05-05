function [T,time,sor_gauss]= steady_sor_gauss(n1,n2,k,T_old,T,dx,dy,X,Y,error,tol,alpha,T_gauss)
    tic                                                                    %start of simulation time
    sor_gauss=1;
    while(error>tol)                                                       %Convergence loop
        for i=2:(n1-1)
            for j=2:(n2-1)
                T(i,j)=((1/(k* dx^2))*(T(i-1,j)+T(i+1,j)))+ ((1/(k* dy^2))*(T(i,j+1)+T(i,j-1)));
                T_gauss=T;
                T(i,j)=((1-alpha)*T_old(i,j)) + (alpha* T_gauss(i,j));
            end
        end
        error=max(max(abs(T-T_old)));                                      %First we are spilitting the matrix into an array by taking the highest values in each of the column and then again max of this array gives the max error value
        sor_gauss=sor_gauss+1;
        T_old=T;
    
        figure(1)   
        [c,d]=contourf(X,Y,T);
        colormap('jet')                                                    %plotting the colormap
        colorbar
        title_text = sprintf('Iteration number=%d for Steady state condition for Successive-Over-Relaxation over Gauss-Seidel Iterative method',sor_gauss);
        title(title_text);
        clabel(c,d)
        xlabel('X')
        ylabel('Y')
        pause(0.03)
    end
    toc;
    time=toc;
end