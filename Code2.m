%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Code 2. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;

% Set up time domain.
Tf = 40.0; T0 = 0.0; Nt = 200; dt = (Tf-T0)/Nt;
tspace = linspace(T0,Tf,Nt+1);

% Set up spatial domain.
Xf = 50.0; X0 = 0.0; Nx = 1000; dx = (Xf-X0)/Nx;
xspace = linspace(X0,Xf,Nx+1);

% Select parameters.
D = 0.5; Estart = -10.0; E0 = -1.0; alpha = 0.5; k0 = 35.0;

% Solve the PDE. (f is a flag for voltammetry type.)
f = 1.0; abspaceL = PDEsolver(xspace,tspace,[D Estart E0 alpha k0 f]);
f = 2.0; abspaceC = PDEsolver(xspace,tspace,[D Estart E0 alpha k0 f]);
f = 3.0; abspaceS = PDEsolver(xspace,tspace,[D Estart E0 alpha k0 f]);

% Post-processing.
IL(1,:) = (abspaceL(2,:)- abspaceL(1,:))/dx; % Get current.
IC(1,:) = (abspaceC(2,:)- abspaceC(1,:))/dx; % Get current.
IS(1,:) = (abspaceS(2,:)- abspaceS(1,:))/dx; % Get current.
EL = VoltE(tspace,Estart,1); % Get voltammetry profile.
EC = VoltE(tspace,Estart,2); % Get voltammetry profile.
ES = VoltE(tspace,Estart,3); % Get voltammetry profile.

%% Figures ...

%%%%%%%%%%%%%%%%%%%%%%%%% Linear sweep voltammetry (DC voltammetry). %%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(xspace,abspaceL(1:Nx+1,round(3.0*Nt/4.0)),'r','linewidth',3.0)
hold on
plot(xspace,abspaceL(1:Nx+1,round(Nt/2.0)),'b','linewidth',3.0)
hold on
plot(xspace,abspaceL(1:Nx+1,round(Nt/4.0)),'Color',[0.0 0.5 0.0],'linewidth',3.0)
hold on
plot(xspace,abspaceL(1:Nx+1,end),'k-','linewidth',3.0)
hold on
plot(xspace,abspaceL(1:Nx+1,1),'k-.','linewidth',3.0)
xlabel('$x$','fontsize',16, 'interpreter','latex')
ylabel('$a(x,t)$','fontsize',16, 'interpreter','latex')
legend({'$t = 10.0$','$t = 20.0$','$t = 30.0$','$t = 40.0$','$t = 0.0$'},...
    'interpreter','latex','Location','best','fontsize',16)
grid on

figure(2)
plot(xspace,abspaceL(Nx+2:2*(Nx+1),round(3.0*Nt/4.0)),'r','linewidth',3.0)
hold on
plot(xspace,abspaceL(Nx+2:2*(Nx+1),round(Nt/2.0)),'b','linewidth',3.0)
hold on
plot(xspace,abspaceL(Nx+2:2*(Nx+1),round(Nt/4.0)),'Color',[0.0 0.5 0.0],'linewidth',3.0)
hold on
plot(xspace,abspaceL(Nx+2:2*(Nx+1),end),'k','linewidth',3.0)
hold on
plot(xspace,abspaceL(Nx+2:2*(Nx+1),1),'k-.','linewidth',3.0)
xlabel('$x$','fontsize',16, 'interpreter','latex')
ylabel('$b(x,t)$','fontsize',16, 'interpreter','latex')
legend({'$t = 10.0$','$t = 20.0$','$t = 30.0$','$t = 40.0$','$t = 0.0$'},...
    'interpreter','latex','Location','best','fontsize',16)
grid on

figure(3)
plot(tspace,IL,'r','linewidth',3.0)
xlabel('$t$','fontsize',16, 'interpreter','latex')
ylabel('$I(t)$','fontsize',16, 'interpreter','latex')
grid on

%%%%%%%%%%%%%%%%%%%%%%%% Cyclic voltammetry. %%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
plot(xspace,abspaceC(1:Nx+1,round(3.0*Nt/4.0)),'r','linewidth',3.0)
hold on
plot(xspace,abspaceC(1:Nx+1,round(Nt/2.0)),'b','linewidth',3.0)
hold on
plot(xspace,abspaceC(1:Nx+1,round(Nt/4.0)),'Color',[0.0 0.5 0.0],'linewidth',3.0)
hold on
plot(xspace,abspaceC(1:Nx+1,end),'k-','linewidth',3.0)
hold on
plot(xspace,abspaceC(1:Nx+1,1),'k-.','linewidth',3.0)
xlabel('$x$','fontsize',16, 'interpreter','latex')
ylabel('$a(x,t)$','fontsize',16, 'interpreter','latex')
legend({'$t = 10.0$','$t = 20.0$','$t = 30.0$','$t = 40.0$','$t = 0.0$'},...
    'interpreter','latex','Location','best','fontsize',16)
grid on

figure(5)
plot(xspace,abspaceC(Nx+2:2*(Nx+1),round(3.0*Nt/4.0)),'r','linewidth',3.0)
hold on
plot(xspace,abspaceC(Nx+2:2*(Nx+1),round(Nt/2.0)),'b','linewidth',3.0)
hold on
plot(xspace,abspaceC(Nx+2:2*(Nx+1),round(Nt/4.0)),'Color',[0.0 0.5 0.0],'linewidth',3.0)
hold on
plot(xspace,abspaceC(Nx+2:2*(Nx+1),end),'k','linewidth',3.0)
hold on
plot(xspace,abspaceC(Nx+2:2*(Nx+1),1),'k-.','linewidth',3.0)
xlabel('$x$','fontsize',16, 'interpreter','latex')
ylabel('$b(x,t)$','fontsize',16, 'interpreter','latex')
legend({'$t = 10.0$','$t = 20.0$','$t = 30.0$','$t = 40.0$','$t = 0.0$'},...
    'interpreter','latex','Location','best','fontsize',16)
grid on

figure(6)
plot(EC,IC,'r','linewidth',3.0)
xlabel('$E(t)$','fontsize',16, 'interpreter','latex')
ylabel('$I(t)$','fontsize',16, 'interpreter','latex')
grid on

%%%%%%%%%%%%%%%%%%%%%%%% Sine wave voltammetry (AC voltammetry). %%%%%%%%%%%%%%%%%%%%%%%%

figure(7)
plot(xspace,abspaceS(1:Nx+1,round(3.0*Nt/4.0)),'r','linewidth',3.0)
hold on
plot(xspace,abspaceS(1:Nx+1,round(Nt/2.0)),'b','linewidth',3.0)
hold on
plot(xspace,abspaceS(1:Nx+1,round(Nt/4.0)),'Color',[0.0 0.5 0.0],'linewidth',3.0)
hold on
plot(xspace,abspaceS(1:Nx+1,end),'k-','linewidth',3.0)
hold on
plot(xspace,abspaceS(1:Nx+1,1),'k-.','linewidth',3.0)
grid on
xlabel('$x$','fontsize',16, 'interpreter','latex')
ylabel('$a(x,t)$','fontsize',16, 'interpreter','latex')
legend({'$t = 10.0$','$t = 20.0$','$t = 30.0$','$t = 40.0$','$t = 0.0$'},...
    'interpreter','latex','Location','best','fontsize',16)
grid on

figure(8)
plot(xspace,abspaceS(Nx+2:2*(Nx+1),round(3.0*Nt/4.0)),'r','linewidth',3.0)
hold on
plot(xspace,abspaceS(Nx+2:2*(Nx+1),round(Nt/2.0)),'b','linewidth',3.0)
hold on
plot(xspace,abspaceS(Nx+2:2*(Nx+1),round(Nt/4.0)),'Color',[0.0 0.5 0.0],'linewidth',3.0)
hold on
plot(xspace,abspaceS(Nx+2:2*(Nx+1),end),'k','linewidth',3.0)
hold on
plot(xspace,abspaceS(Nx+2:2*(Nx+1),1),'k-.','linewidth',3.0)
grid on
xlabel('$x$','fontsize',16, 'interpreter','latex')
ylabel('$b(x,t)$','fontsize',16, 'interpreter','latex')
legend({'$t = 10.0$','$t = 20.0$','$t = 30.0$','$t = 40.0$','$t = 0.0$'},...
    'interpreter','latex','Location','best','fontsize',16)
grid on

figure(9)
plot(EL,IS,'r','linewidth',3.0)
xlabel('$E_{DC}(t)$','fontsize',16, 'interpreter','latex')
ylabel('$I(t)$','fontsize',16, 'interpreter','latex')
grid on

%% Functions ...

% Voltammetry functions.

function [E] = VoltE(t,Estart,f)

    if (f == 1) % Linear voltammetry.
    
        E = t + Estart;
    
    elseif (f == 2) % Cyclic voltammetry.
        
        trev = 20.0;
        v = mod(floor(t/trev),2);
        E(v==0) = Estart + t(v==0) - floor(t(v==0)/trev)*trev;
        E(v==1) = Estart + 2.0*trev - t(v==1) + (floor(t(v==1)/trev)-1.0)*trev;
        
    elseif (f == 3) % Sine voltammetry.
        
        omega = 16.0; DE = 2.0;
        E = Estart + t + DE*sin(omega*t);
        
    else
        
       fprintf('Unsupported voltammetry type.\n')
       
    end
    
end

function [abspace] = PDEsolver(xspace,tspace,par)

    % Open parameter vector.
    D = par(1,1); Estart = par(1,2); E0 = par(1,3); alpha = par(1,4); k0 = par(1,5); f = par(1,6);
    
    % Create temporal and spatial grid.    
    Nt = size(tspace,2)-1; dt = tspace(1,2)-tspace(1,1);
    Nx = size(xspace,2)-1; dx = xspace(1,2)-xspace(1,1);

    % Set up ab-space.
    abspace = zeros(2*(Nx+1),Nt+1);
    % Implement ICs.
    abspace(1:Nx+1,1) = 1.0; abspace(Nx+2:2*(Nx+1),1) = 0.0; % Initial conditions.

    % Select diffusion constant and theta.
    theta = 0.5; mu = dt/(dx)^2;
    
    % Set up the matrices.
    e = ones(Nx+1,1);
    A1 = spdiags([(1.0-theta)*mu*e (1.0-2.0*(1.0-theta)*mu)*e (1.0-theta)*mu*e],-1:1,Nx+1,Nx+1);
    B1 = spdiags([(1.0-theta)*D*mu*e (1.0-2.0*(1.0-theta)*D*mu)*e (1.0-theta)*D*mu*e],-1:1,Nx+1,Nx+1);
    A2 = spdiags([-theta*mu*e (1.0+2.0*theta*mu)*e -theta*mu*e],-1:1,Nx+1,Nx+1);
    B2 = spdiags([-theta*mu*D*e (1.0+2.0*theta*D*mu)*e -theta*mu*D*e],-1:1,Nx+1,Nx+1);

    % Implement time independent BCs.
    M1 = blkdiag(A1,B1); M2 = blkdiag(A2,B2);
    M1(1,:) = 0.0; M1(Nx+1,:) = 0.0; M1(Nx+2,:) = 0.0; M1(2*(Nx+1),:) = 0.0;
    M2(1,:) = 0.0; M2(Nx+1,:) = 0.0; M2(Nx+2,:) = 0.0; M2(2*(Nx+1),:) = 0.0;
    M2(Nx+2,1) = -1.0; M2(Nx+2,2) = 1.0;  M2(Nx+2,Nx+2) = -1.0; M2(Nx+2,Nx+3) = 1.0; % Conservation of electrons.
    M2(Nx+1,Nx+1) = 1.0; M2(2*(Nx+1),2*(Nx+1)) = 1.0; % Far field conditions.
    m = zeros(2*(Nx+1),1); m(Nx+1,1) = 1.0;

    % March in time.
    for i = 1:Nt
        % Implement time dependent BCs.
        M2(1,1) = -1.0/dx-k0*exp((1.0-alpha)*(VoltE(tspace(1,i),Estart,f)-E0));
        M2(1,2) = 1.0/dx;
        M2(1,Nx+2) = k0*exp(-alpha*(VoltE(tspace(1,i),Estart,f)-E0));

        abspace(:,i+1) = M2\(M1*abspace(:,i)+m);    
    end

end

