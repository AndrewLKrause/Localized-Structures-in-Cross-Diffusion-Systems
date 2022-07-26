clear; close all;


a=0.075; b=0.029;
d1 = 1;  d4=0.9; d3 = 0.6; d2 = 1;

%Domain length (needs to be reasonably large).
L = 50;

%Number of grid points.
N = 1000;

dx = L/N;

%Kinetic functions.
f = @(u,v)(a-u+u.^2.*v);
g = @(u,v)(b - u.^2.*v);

%Steady state values.
%uss = a+b; vss=b/(a+b)^2;

%
%

%localized init condition
uss = a+b; vss = b/uss^2;
x = linspace(0,L,N);
%initial conditions 

spike = 3*exp(-1*x.^2);
uss=uss+spike;
vss=vss+spike;

%uss = a+b+1e-2*randn(N,1);vss=b./((a+b).^2) +1e-2*randn(N,1);

uinit = [uss;vss];

%Create the negative Laplacian
%[~,~,lap] = laplacian([N],{'NN'});
%lap = -lap*(N/L)^2;
lap = spdiags([ones(N,1),-2*ones(N,1),ones(N,1)],[1,0,-1],N,N);
%Neumann conditions:
lap(1) = -1; lap(end) = -1;
lap = (1/dx)^2*lap;



%Matlab ODE solvers need 
F = @(t,U)[f(U(1:N),U(N+1:2*N))+d1*(1/dx)^2*Lap(U(1:N))+d2*(1/dx)^2*Lap(U(N+1:2*N));g(U(1:N),U(N+1:2*N))+d3*(1/dx)^2*Lap(U(1:N))+d4*(1/dx)^2*Lap(U(N+1:2*N))];

%Get solution at times
T = linspace(0,15000,20000);

%Create Jacobian sparsity patern.
LP = lap~=0; I = eye(N);
JPattern = [LP, I+LP; I+LP, LP];

%Solve the system using a stiff solver and low tolerances.
odeparams = odeset('RelTol',1e-6,'AbsTol',1e-6,'JPattern',JPattern);
[T, U] = ode15s(F,T,uinit,odeparams);

u = U(:,1:N); v = U(:,N+1:2*N);
%close all;
%imagesc(u);colorbar
plot(u(end,:))


function u = Lap(u)
    %u(i-1)-2*u(i)+u(i+1)
    u = [u(2)-u(1); u(1:end-2)-2*u(2:end-1)+u(3:end);u(end-1)-u(end)];
end