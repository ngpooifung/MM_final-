%% Code 3 %%
%clear all; clc;

% Set up time domain.
Tf = 40.0; T0 = 0.0; Nt = 200; dt = (Tf-T0)/Nt;
tspace = linspace(T0,Tf,Nt+1);

% Set up spatial domain.
Xf = 50.0; X0 = 0.0; Nx = 1000; dx = (Xf-X0)/Nx;
xspace = linspace(X0,Xf,Nx+1);

% Select intensity of noise
r = 0.05;

% Number of iterations
N = 5;

% Select control parameters.
D = 0.5; Estart = -10.0;

% Select parameters to be inferred
E0 = 0; alpha = 0.5; k0 = 0.1;

% Create experimental data. (f is a flag for voltammetry type.)
f = 1.0; abspaceL = PDEsolver(xspace,tspace,[D Estart E0 alpha k0 f]);
f = 2.0; abspaceC = PDEsolver(xspace,tspace,[D Estart E0 alpha k0 f]);
f = 3.0; abspaceS = PDEsolver(xspace,tspace,[D Estart E0 alpha k0 f]);
% Get current from experimental data
IL(1,:) = (abspaceL(2,:)- abspaceL(1,:))/dx; 
IC(1,:) = (abspaceC(2,:)- abspaceC(1,:))/dx; 
IS(1,:) = (abspaceS(2,:)- abspaceS(1,:))/dx; 

% Iterations to retreive paramets
for ii = 1:N
    % Add noise to data
    IL(1,:) = IL(1,:) + r*max(IL)*randn(1,length(IL)); 
    IC(1,:) = IC(1,:) + r*max(IC)*randn(1,length(IC)); 
    IS(1,:) = IS(1,:) + r*max(IS)*randn(1,length(IS));
    
    % Minimise error with contrained optimisation
    f = 1; paramL = fminsearch(@(par1) minfun(par1,[D Estart f dx],xspace,...
        tspace,IL), [-1 0.2 2]);
    f = 2; paramC = fminsearch(@(par1) minfun(par1,[D Estart f dx],xspace,...
        tspace,IC), [-1 0.2 2]);
    f = 3; paramS = fminsearch(@(par1) minfun(par1,[D Estart f dx],xspace,...
        tspace,IS), [-1 0.2 2]);
    
    % Save results
    paramsL(ii, :) = paramL;
    paramsC(ii, :) = paramC;
    paramsS(ii, :) = paramS;
end
colNames = {(['E0 = ',num2str(E0)]),(['alpha = ',num2str(alpha)]),(['k0 = ',num2str(k0)])};
paramsLTab = array2table(paramsL,'VariableNames',colNames)
paramsCTab = array2table(paramsC,'VariableNames',colNames)
paramsSTab = array2table(paramsS,'VariableNames',colNames)

params = [E0 alpha k0];
for jj = 1:3
    errL = norm(paramsL(:,jj) - params(jj));
    errC = norm(paramsC(:,jj) - params(jj));
    errS = norm(paramsS(:,jj) - params(jj));
    
    ErrL(jj) = errL; ErrC(jj) = errC; ErrS(jj) = errS; 
end
colNames = {'error E0','error alpha','error k0'};
ErrLTab = array2table(ErrL,'VariableNames',colNames)
ErrCTab = array2table(ErrC,'VariableNames',colNames)
ErrSTab = array2table(ErrS,'VariableNames',colNames)

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

function [errsqrd] = minfun(par1,par2,xspace,tspace,expdataI)
% Function to be minimised - squared error between experimental data
% (input) and numerical data

% Get parameters
E0 = par1(1,1); alpha = par1(1,2); k0 = par1(1,3);
D = par2(1,1); Estart = par2(1,2); f = par2(1,3); dx = par2(1,4);

% Calculate error
if alpha <= 0 || alpha >= 1 || k0 < 0 % If parameters move outside allowed 
    errsqrd = 100000;                % region, make error large
else
    abspace = PDEsolver(xspace,tspace,[D Estart E0 alpha k0 f]);
    I(1,:) = (abspace(2,:)- abspace(1,:))/dx;
    errsqrd = norm(expdataI - I);
end
end
