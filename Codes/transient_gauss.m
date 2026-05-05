function [T,time,gauss_iter]= transient_gauss(n1,n2,T_old,T_initial,T,dx,dy,x,y,tol,term1,term2,term3,gauss_iter)
tic                                                                        %start of simulation time


for p=1:1400                                                               %time loop
    error=9e9;
    while(error>tol)                                                       %convergence loop
    
        
            for i=2:(n1-1)
                for j=2:(n2-1)
                    H=(T(i-1,j) + T_old(i+1,j));
                    V=(T(i,j-1) + T_old(i,j+1));
                    T(i,j)= (term1*T_initial(i,j)) + (term2*H) + (term3*V);
                    
                end
            end
       
            error=max(max(abs(T-T_old)));                                  %First we are spilitting the matrix into an array by taking the highest values in each of the column and then again max of this array gives the max error value
            T_old=T;                                                       %Upddating the values after each iteration
            gauss_iter=gauss_iter+1;    
    end
            T_initial = T;                                                 %updating the value after each time step
            figure(1)   
            [c,d]=contourf(x,y,T);
            colormap('jet');                                               %plotting the colormap
            colorbar;
            clabel(c,d)
            title_text = sprintf('Iteration number=%d for Steady state condition using Gauss-Seidel Implicit Iterative Solver',gauss_iter);
            title(title_text);
            xlabel('X')
            ylabel('Y')
            pause(0.03)
    
    end
    toc;                                                                   %end time of simulation
    time=toc;
end