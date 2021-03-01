% Solution of 1D diffusion-convection problem for Dirichlet BC's using the
% linear approximation for diffusion and upwind for the convection

clear all;
close all;

% Inputs
rho=1;            
uinf=1;    
L=1;              
N=365;            
gamma=0.1;
dx=L/N;           
u1=1;
un=0;
flux=rho*uinf;

tic
for j=2:1:N-1       %(Loop on number of cells that are not on boundary)
                    %(Since global axes in this example will contain
                    %a_{i-1}, a_{i+1} we can not use the loop outside 2:N-1
    ae=gamma/dx;        %Value of ae in local axes
    aw=gamma/dx+flux;        %Value of aw in local axes
    ap=aw+ae;       %Value of ap in local axes
    
    a(j,j-1)=-aw;   %Values in global axes
    a(j,j)=ap;      %Values in global axes
    a(j,j+1)=-ae;   %Values in global axes
    b(j,1)=0;       %Value in source term
end

% In this code the boundaries are hardwired.
j=1;
aw=2*gamma/dx+flux;
ae=gamma/dx;
ap=flux+gamma/dx+2*gamma/dx;
       
a(j,j)=ap;
a(j,j+1)=-ae;
b(j,1)=aw*u1;

% In this code the boundaries are hardwired.
j=N;
aw=gamma/dx+flux;
ae=2*gamma/dx;
ap=2*gamma/dx+gamma/dx+flux;

a(j,j-1)=-aw;
a(j,j)=ap;
b(j,1)=ae*un;

% End Assembly

sol=a\b;
xsol=dx/2:dx:L-dx/2;
figure;
plot(xsol,sol,'r')
xlabel('axial location x');ylabel('\phi-distribution')

%The exact solution
xexact=xsol;
temp=(exp(flux*xexact/gamma)-1)./(exp(flux*L/gamma)-1);
uexact=u1+(un-u1)*temp;
hold;plot(xexact,uexact,'b')
legend('Numerical','Exact')

format long
errorUpw=max(abs(sol'-uexact))
toc    

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
