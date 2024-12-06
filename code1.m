Xf = 50.0; X0 = 0.0; Nx = 1000; dx = (Xf-X0)/Nx;
xspace = linspace(X0,Xf,Nx+1);
Tf = 40.0; T0 = 0.0; Nt = 200; dt = dx^2/2;
tspace = linspace(T0,Tf,Nt+1);
mu=dt/(dx^2);
theta=1/2;
t=[10,20,30,40,50];

e=ones(Nx+1,1);
A=spdiags([e -2*e e],-1:1,Nx+1,Nx+1);
A(1,:)=0;
A(end,:)=0;
I=sparse(eye(Nx+1));
I2=I;
I2(1,1)=0;
I2(end,end)=0;
g=zeros(Nx+1,1);
g(end)=1;
P=I-mu*theta*A;
Q=I2+mu*(1-theta)*A;
a2=zeros(Nx+1,length(t));
b2=zeros(Nx+1,length(t));
c2=zeros(Nx+1,length(t));
% a3=zeros(Nx+1,Nt+1);
% a=ones(Nx+1,1);
% a3(:,1)=a;
% 
% for i=1:Nt
%     a=P\(Q*a+g);
%     a3(:,i+1)=a;
% end
% 
% C=(a3(2,:)-a3(1,:))/dx;

for i=1:length(t)
    a=ones(Nx+1,1);
    for j=1:(t(i)/dt)
        a=P\(Q*a+g);
    end
    a2(:,i)=a;
    b2(:,i)=ones(Nx+1,1)-a;
end

for i=1:length(t)
    c2(:,i)=erf(xspace./(2*sqrt(t(i))));
end


figure
plot(xspace,ones(Nx+1,1),'--k','Linewidth',3)
hold on 
plot(xspace,a2(:,1),xspace,a2(:,2),xspace,a2(:,3),xspace,a2(:,4),xspace,a2(:,5),'LineWidth',3)
legend({'$t = 0$','$t = 10$','$t = 20$','$t = 30$','$t = 40$','$t = 50$',},...
    'interpreter','latex','Location','best','fontsize',16)
grid on
xlim([0 50])
ylim([0 1.2])
xlabel('$x$','fontsize',16, 'interpreter','latex')
ylabel('$a(x,t)$','fontsize',16, 'interpreter','latex')

figure
plot(xspace,zeros(Nx+1,1),'--k','Linewidth',3)
hold on 
plot(xspace,b2(:,1),xspace,b2(:,2),xspace,b2(:,3),xspace,b2(:,4),xspace,b2(:,5),'LineWidth',3)
legend({'$t = 0$','$t = 10$','$t = 20$','$t = 30$','$t = 40$','$t = 50$'},...
    'interpreter','latex','Location','best','fontsize',16)
grid on
xlim([0 50])
ylim([0 1.2])
xlabel('$x$','fontsize',16, 'interpreter','latex')
ylabel('$b(x,t)$','fontsize',16, 'interpreter','latex')

figure
plot(xspace,abs(c2(:,1)-a2(:,1)),xspace,abs(c2(:,2)-a2(:,2)),xspace,abs(c2(:,3)-a2(:,3)),xspace,abs(c2(:,4)-a2(:,4)),xspace,abs(c2(:,5)-a2(:,5)),'LineWidth',3)
legend({'$t = 10$','$t = 20$','$t = 30$','$t = 40$','$t = 50$'},...
    'interpreter','latex','Location','best','fontsize',16)
xlabel('$x$','fontsize',16, 'interpreter','latex')
ylabel('error in $a$','fontsize',16, 'interpreter','latex')
grid on

% figure
% plot(tspace,C,'LineWidth',3)
% xlabel('$t$','fontsize',16, 'interpreter','latex')
% ylabel('Current','fontsize',16, 'interpreter','latex')
% grid on

