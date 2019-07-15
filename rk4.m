clc
%RK4
%declare variable
a=0;
b=4;
h=0.001;
N=(b-a)/h;
%higher order to ODE
g1=@(eta,f,p,q) p;
g2=@(eta,f,p,q) q;
g3=@(eta,f,p,q) -f*q+2*p^(2);
%initial value problem
eta(1)=a;
f(1)=0;
p(1)=1;
%initial guess
q(1)=-1.284367;

for i=1:N
    k1=g1(eta(i),f(i),p(i),q(i));
    m1=g2(eta(i),f(i),p(i),q(i));
    n1=g3(eta(i),f(i),p(i),q(i));
    
    k2=g1(eta(i)+h/2,f(i)+h*k1/2,p(i)+h*m1/2,q(i)+h*n1/2);
    m2=g2(eta(i)+h/2,f(i)+h*k1/2,p(i)+h*m1/2,q(i)+h*n1/2);
    n2=g3(eta(i)+h/2,f(i)+h*k1/2,p(i)+h*m1/2,q(i)+h*n1/2);
    
    k3=g1(eta(i)+h/2,f(i)+h*k2/2,p(i)+h*m2/2,q(i)+h*n2/2);
    m3=g2(eta(i)+h/2,f(i)+h*k2/2,p(i)+h*m2/2,q(i)+h*n2/2);
    n3=g3(eta(i)+h/2,f(i)+h*k2/2,p(i)+h*m2/2,q(i)+h*n2/2);
    
    k4=g1(eta(i)+h,f(i)+h*k3,p(i)+h*m3,q(i)+h*n3);
    m4=g2(eta(i)+h,f(i)+h*k3,p(i)+h*m3,q(i)+h*n3);
    n4=g3(eta(i)+h,f(i)+h*k3,p(i)+h*m3,q(i)+h*n3);
     
    f(i+1)=f(i)+h*(k1+2*k2+2*k3+k4)/6;
    p(i+1)=p(i)+h*(m1+2*m2+2*m3+m4)/6;
    q(i+1)=q(i)+h*(n1+2*n2+2*n3+n4)/6;
    eta(i+1)=eta(i)+h;
    
    %print output
    fprintf('Iteration=%2d\t eta(%d)=%.4f\t f(%d)=%.4f\t p(%d)=%.4f\t q(%d)=%.4f\t \n',i,i,eta(i+1),i,f(i+1),i,p(i+1),i,q(i));
%display output in excel
     xlswrite('A1RK40.01-result.xlsx',[eta.',f.',p.',q.']);
end

%plot graph
%since x-axis is actually an eta, so we need to plot the graph f,p,q with
%respect to eta
plot(eta,f, 'LineWidth', 1)
hold on
plot(eta,p, 'LineWidth', 1)
plot(eta,q, 'LineWidth', 1)
axis([0 4 -2 inf])
hold off
legend('f', 'p', 'q')
title('Runge Kutta') 
xlabel('\eta')
ylabel('f,p,q')


    


