 clc
%declare variable
a=0;
b=4;
h=0.01;
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
%RKF
for i=1:N
    k1=g1(eta(i),f(i),p(i),q(i));
    m1=g2(eta(i),f(i),p(i),q(i));
    n1=g3(eta(i),f(i),p(i),q(i));
   
    k2=g1(eta(i)+h/4,f(i)+h*k1/4,p(i)+h*m1/4,q(i)+h*n1/4);
    m2=g2(eta(i)+h/4,f(i)+h*k1/4,p(i)+h*m1/4,q(i)+h*n1/4);
    n2=g3(eta(i)+h/4,f(i)+h*k1/4,p(i)+h*m1/4,q(i)+h*n1/4);
    
    k3=g1(eta(i)+3*h/8,f(i)+3*h*k1/32+9*h*k2/32,p(i)+3*h*m1/32+9*h*m2/32,q(i)+3*h*n1/32+9*h*n2/32);
    m3=g2(eta(i)+3*h/8,f(i)+3*h*k1/32+9*h*k2/32,p(i)+3*h*m1/32+9*h*m2/32,q(i)+3*h*n1/32+9*h*n2/32);
    n3=g3(eta(i)+3*h/8,f(i)+3*h*k1/32+9*h*k2/32,p(i)+3*h*m1/32+9*h*m2/32,q(i)+3*h*n1/32+9*h*n2/32);
    
    k4=g1(eta(i)+12*h/13,f(i)+1932*h*k1/2197-7200*h*k2/2197+7296*h*k3/2197,p(i)+1932*h*m1/2197-7200*h*m2/2197+7296*h*m3/2197,q(i)+1932*h*n1/2197-7200*h*n2/2197+7296*h*n3/2197);
    m4=g2(eta(i)+12*h/13,f(i)+1932*h*k1/2197-7200*h*k2/2197+7296*h*k3/2197,p(i)+1932*h*m1/2197-7200*h*m2/2197+7296*h*m3/2197,q(i)+1932*h*n1/2197-7200*h*n2/2197+7296*h*n3/2197);
    n4=g3(eta(i)+12*h/13,f(i)+1932*h*k1/2197-7200*h*k2/2197+7296*h*k3/2197,p(i)+1932*h*m1/2197-7200*h*m2/2197+7296*h*m3/2197,q(i)+1932*h*n1/2197-7200*h*n2/2197+7296*h*n3/2197);
    
    k5=g1(eta(i)+h,f(i)+439*h*k1/216-8*h*k2+3680*h*k3/513-845*h*k4/4104,p(i)+439*h*m1/216-8*h*m2+3680*h*m3/513-845*h*m4/4104,q(i)+439*h*n1/216-8*h*n2+3680*h*n3/513-845*h*n4/4104);
    m5=g2(eta(i)+h,f(i)+439*h*k1/216-8*h*k2+3680*h*k3/513-845*h*k4/4104,p(i)+439*h*m1/216-8*h*m2+3680*h*m3/513-845*h*m4/4104,q(i)+439*h*n1/216-8*h*n2+3680*h*n3/513-845*h*n4/4104);
    n5=g3(eta(i)+h,f(i)+439*h*k1/216-8*h*k2+3680*h*k3/513-845*h*k4/4104,p(i)+439*h*m1/216-8*h*m2+3680*h*m3/513-845*h*m4/4104,q(i)+439*h*n1/216-8*h*n2+3680*h*n3/513-845*h*n4/4104);
   
    k6=g1(eta(i)+h/2,f(i)-8*h*k1/27+2*h*k2-3554*h*k3/2565+1859*h*k4/4104-11*h*k5/40,p(i)-8*h*m1/27+2*h*m2-3554*h*m3/2565+1859*h*m4/4104-11*h*m5/40,q(i)-8*h*n1/27+2*h*n2-3554*h*n3/2565+1859*h*n4/4104-11*h*n5/40);
    m6=g2(eta(i)+h/2,f(i)-8*h*k1/27+2*h*k2-3554*h*k3/2565+1859*h*k4/4104-11*h*k5/40,p(i)-8*h*m1/27+2*h*m2-3554*h*m3/2565+1859*h*m4/4104-11*h*m5/40,q(i)-8*h*n1/27+2*h*n2-3554*h*n3/2565+1859*h*n4/4104-11*h*n5/40);
    n6=g3(eta(i)+h/2,f(i)-8*h*k1/27+2*h*k2-3554*h*k3/2565+1859*h*k4/4104-11*h*k5/40,p(i)-8*h*m1/27+2*h*m2-3554*h*m3/2565+1859*h*m4/4104-11*h*m5/40,q(i)-8*h*n1/27+2*h*n2-3554*h*n3/2565+1859*h*n4/4104-11*h*n5/40);
    
    f(i+1)=f(i)+h*(125*k1/216+7040*k3/2565+10985*k4/4101-k5)/5;
    p(i+1)=p(i)+h*(125*m1/216+7040*m3/2565+10985*m4/4101-m5)/5;
    q(i+1)=q(i)+h*(125*n1/216+7040*n3/2565+10985*n4/4101-n5)/5;
    eta(i+1)=eta(i)+h;
    
    %print output
    fprintf('Iteration=%2d\t eta(%d)=%.4f\t f(%d)=%.4f\t p(%d)=%.4f\t q(%d)=%.4f\t \n',i,i,eta(i+1),i,f(i+1),i,p(i+1),i,q(i));
   %display output in excel
 xlswrite('A1RKf0.01-result.xlsx',[eta.',f.',p.',q.']);
end

%plot graph
%since x-axis is actpally an eta, so we need to plot the graph f,p,q with
%respect to eta
plot(eta,f, 'LineWidth', 1)
hold on
plot(eta,p, 'LineWidth', 1)
plot(eta,q, 'LineWidth', 1)
axis([0 4 -2 inf])
hold off
legend('f', 'p', 'q')
title('Runge Kutta Fehlberg') 
xlabel('\eta')
ylabel('f,p,q')


