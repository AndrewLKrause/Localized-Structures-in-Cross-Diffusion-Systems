clear; close all;


d1=1; d4=1; rho=1;chi0=4.01;

%Domain length (needs to be reasonably large).
L = 50;

%Number of grid points.
N = 1000;

dx = L/N;

%Kinetic functions.
f = @(u)(u-u.^2);
g = @(u,v)(rho*u-v);

%Homogeneous background steady state
uss = 1; vss = rho*uss;
x = linspace(0,L,N);


%initial conditions 
rng(1);
uss=uss*(1+1e-1*randn(N,1));
vss=vss*(1+1e-1*randn(N,1));

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
F = @(t,U)[f(U(1:N))+d1*(1/dx)^2*Lap(U(1:N))-chi0*(1/(2*dx^2))*chmtxis(U(1:N),U(N+1:2*N));...,
    g(U(1:N),U(N+1:2*N))+d4*(1/dx)^2*Lap(U(N+1:2*N))];

%Get solution at times
T = linspace(0,15000,20000);

%Create Jacobian sparsity patern.
LP = lap~=0; I = eye(N);
JPattern = [LP, I+LP; I+LP, LP];

%Solve the system using a stiff solver and low tolerances.
odeparams = odeset('RelTol',1e-11,'AbsTol',1e-11,'JPattern',JPattern);
[T, U] = ode15s(F,T,uinit,odeparams);

u = U(:,1:N); v = U(:,N+1:2*N);
%close all;
%imagesc(u);colorbar
plot(u(end,:))


function u = Lap(u)
    u = [u(2)-u(1); u(1:end-2)-2*u(2:end-1)+u(3:end);u(end-1)-u(end)];
end


%Version of (u_(i+1)+u_i)*(v_(i+1)-v_i)-(u_i+u_(i-1))*(v_i-v_(i-1))
function chmtxis = chmtxis(u,v)
    chmtxis = [(u(2)+u(1)).*(v(2)-v(1));...,
        (u(3:end)+u(2:end-1)).*(v(3:end)-v(2:end-1))-(u(2:end-1)+u(1:end-2)).*(v(2:end-1)-v(1:end-2));...,
        -(u(end)+u(end-1)).*(v(end)-v(end-1))];
end