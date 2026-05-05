function [T,time,sor_jac]= steady_sor_jac(n1,n2,k,T_old,T,dx,dy,X,Y,error,tol,alpha,T_jac)
tic
sor_jac=1;
    while(error>tol)
        for i=2:(n1-1)
            for j=2:(n2-1)
                T(i,j)=((1/(k* dx^2))*(T_old(i-1,j)+T_old(i+1,j)))+ ((1/(k* dy^2))*(T_old(i,j+1)+T_old(i,j-1)));
                T_jac=T;
                T(i,j)=((1-alpha)*T_old(i,j)) + (alpha* T_jac(i,j));
            end
        end
        error=max(max(abs(T-T_old)));
        sor_jac=sor_jac+1;
        T_old=T;
        figure(1)   
        [c,d]=contourf(X,Y,T);
        colormap('jet')
        colorbar
        title_text = sprintf('Iteration number=%d for Steady state condition for Successive-Over-Relaxation over Jacobian Iterative method',sor_jac);
        title(title_text);
        clabel(c,d)
        xlabel('X')
        ylabel('Y')
        pause(0.03)
    end
    toc;
    time=toc;
end