function [T,time,exp_iter]=transient__explicit_1(n1,n2,T_old,T_initial,T,X,Y,tol,k1,k2,exp_iter)
tic                                                                        %start of simulation time

    for p=1:1400                                                           %time loop
        error=9e9;
        while(error>tol)                                                   %convergence loop
    
        
            for i=2:(n1-1)
                for j=2:(n2-1)
                    H=(T_old(i-1,j) - 2*T_old(i,j) + T_old(i+1,j));
                    V=(T_old(i,j-1) - 2*T_old(i,j) + T_old(i,j+1));
                    T(i,j)= T_initial(i,j) + (k1*H) + (k2*V);
                    
                end
            end
       
            error=max(max(abs(T-T_old)));                                  %First we are spilitting the matrix into an array by taking the highest values in each of the column and then again max of this array gives the max error value
            T_old=T;                                                       %Upddating the values after each iteration
            exp_iter=exp_iter+1;
            
         end
            T_initial = T;                                                 %updating the value after each time step
            figure(1)   
            [c,d]=contourf(X,Y,T);
            colormap('jet');                                               %plotting the colormap
            colorbar;
            clabel(c,d);
            title_text = sprintf('Iteration number=%d for Transient state condition using Explicit Iterative Solver',exp_iter);
            title(title_text);
            xlabel('X')
            ylabel('Y')
            pause(0.03)
    
    end
    toc;                                                                   %end time of simulation
    time=toc;
end 