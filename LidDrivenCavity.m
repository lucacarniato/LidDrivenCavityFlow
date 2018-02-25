% Lid driven cavity: calculates the steady state flow in a cavity where the top wall slides
% plots the velocity fields and the reference solution by ghia et al. 

clc;
clear;
close all;

Re       = 100;  %Reynolds number
N        = 12;   %Number of volumes in x and y directions
dt       = 0.07; %timestep
Delta    = 1/N;  %disc 
maxIter  = 5000; %max number explicit iterations 
divTol   = 1e-5; %exit if flow divergence is larger than this number
ssTol    = 1e-3; %steady state criteria
rho      = 1;    %density
relaxFac = 0.4;
upwind   = 1;
implicit = 1; 
simple   = 1;   
debug    = 0;

%mesh generation
x = zeros(N+1,1);
for i=1:N+1
    xi = (i-1)*Delta;
    %non-uniform
    %x(i) = 0.5*(1. - cos(pi*xi)); 
    x(i) = (i-1)*Delta;
end

h = zeros(N+2,1);
for i=1:N
    h(i+1) = x(i+1)-x(i);
end

%tangential velocites
for i=1:N+2
    uUP(i) =  1; %punctual velocity 
    uLO(i) =  0;
    vLF(i) =  0; 
    vRH(i) =  0;
end

%set up pressure matrix
OFFSET  = 1;
for i=1:N
    for j=1:N
        k = (i-1)*N + j;
        
        APR=-2*h(i+OFFSET)*dt/((h(j+OFFSET  )+h(j+1+OFFSET))*rho); %u i + 0.5
        BPR=-2*h(i+OFFSET)*dt/((h(j-1+OFFSET)+h(j+OFFSET  ))*rho); %u i - 0.5
        CPR=-2*h(j+OFFSET)*dt/((h(i+OFFSET  )+h(i+1+OFFSET))*rho); %v i + 0.5
        DPR=-2*h(j+OFFSET)*dt/((h(i-1+OFFSET)+h(i+OFFSET  ))*rho); %v i - 0.5
      
        if i == 1
            if j==1
                B(k,3) = -APR-CPR;
                B(k+1,4) = APR;
                B(k+N,5) = CPR;
            elseif j==N
                B(k-1,2) = BPR;
                B(k,3)   = -BPR-CPR;
                B(k+N,5) = CPR;
            else
                B(k-1,2) = BPR;
                B(k,3)   = -APR-BPR-CPR;
                B(k+1,4) = APR;
                B(k+N,5) = CPR;
            end;
        elseif i==N
            if j==1
                B(k-N,1) = DPR;
                B(k,3) = -APR-DPR;
                B(k+1,4) = APR;
            elseif j==N
                B(k-N,1) = DPR;
                B(k-1,2) = BPR;
                B(k,3) = -BPR-DPR;
            else
                B(k-N,1) = DPR;
                B(k-1,2) = BPR;
                B(k,3) = -APR-BPR-DPR;
                B(k+1,4) = APR;
            end
        elseif j == 1
            if ( i > 1 & i < N )
                B(k-N,1) = DPR;
                B(k,3)   = -APR-CPR-DPR;
                B(k+1,4) = APR;
                B(k+N,5) = CPR;
            end
        elseif j == N
            if ( i > 1 & i < N )
                B(k-N,1) = DPR;
                B(k-1,2) = BPR;
                B(k,3)   = -BPR-CPR-DPR;
                B(k+N,5) = CPR;
            end
        else
            B(k-N,1) = DPR;
            B(k-1,2) = BPR;
            B(k,3)   = -APR-BPR-CPR-DPR;
            B(k+1,4) = APR;
            B(k+N,5) = CPR;
        end
    end
end
            
d = [ -N
    -1
    0
    1
    N];

A = spdiags(B,d,N*N,N*N);

%Initialization
iter=0;
udiff=1;
vdiff=1;

pres  = zeros(N,N);
pvec  = zeros(N*N,1);
pc    = zeros(N,N);
pcvec = zeros(N*N,1);
diver = zeros(N,N);
f     = zeros(N*N,1);
Ru    = zeros(N,N+1);
Rv    = zeros(N+1,N);

%hold coefficents for u and v
diagConst      = zeros((N+2)*(N+1),5);
diagConst(:,3) = ones((N+2)*(N+1),1);
diag           = diagConst;
diagPos        = [-(N+2) -1 0 1 (N+2)];
rhsConst       = zeros((N+2)*(N+1),1);
rhs            = rhsConst;
vsol           = zeros((N+2)*(N+1),1);   
vexp           = zeros((N+2)*(N+1),1);  
np_upwind      = zeros(N,1);
n_upwind       = zeros(N,1);
tp_upwind      = zeros(N,1);
t_upwind       = zeros(N,1);
np_upwind(:) =0.5; 
n_upwind(:)  = 0.5;
tp_upwind(:) = 0.5;
t_upwind(:)  = 0.5;
%normal
np_interface = zeros(N,1);
n_interface  = zeros(N,1);
%tangential
tp_interface = zeros(N,1);
t_interface  = zeros(N,1);

u  = zeros(N+2,N+1);
v  = zeros(N+1,N+2);
initval = 0;
if (debug == 1)
u  = ones(N+2,N+1);
v  = ones(N+1,N+2);
initval =1;
end
uold = u; 
vold = v;

% grid terms that do not change
adv          = zeros(N,3);
normalstress = zeros(N,2);
tgstress     = zeros(N,4);
tgstress(:,1) =[1./(Re*h(3:N+1).* (h(2:N)+h(3:N+1)));      0];
tgstress(:,2) =[1./(Re*h(2:N  ).* (h(2:N)+h(3:N+1)));      0];
tgstress(:,3) =[0      ;1./(Re*h(3:N+1).* (h(3:N+1)+h(2:N)))];
tgstress(:,4) =[0      ;1./(Re*h(2:N).* (h(3:N+1)+h(2:N)))  ];
dimedis       = (rho* (h(2:N)+h(3: N+1)))./(2 * dt);
if (~implicit)
dimedis           = -dimedis;
end

for j = 2: N
  adv(j,:) = 0.25 * rho  *[h(j-1)+h(j); h(j)+h(j+1);h(j+1)+h(j+2)];
  normalstress(j,:) = [2/(Re*h(j+1));2/(Re*h(j))];
  diagInd = ((j-1)*(N + 2) + 2) : ((j-1)*(N + 2) + (N+1));
  %set the members that do not change
  diagConst(diagInd-(N+2),1) = - (h(j)+h(j+1)).* tgstress(:,4);
  diagConst(diagInd-1,2)     = - normalstress(j,1);
  diagConst(diagInd,3)       = dimedis(j-1) + normalstress(j,1)+normalstress(j,2) ...
  + (h(j)+h(j+1)).*(tgstress(:,2)+tgstress(:,3));
  diagConst(diagInd+1,4)     =  - normalstress(j,2);
  diagConst(diagInd+(N+2),5) =  - (h(j)+h(j+1)).* tgstress(:,1);
  diagConst(diagInd(1),  3) = diagConst(diagInd(1),  3) + (h(j)+h(j+1))/(Re*(h(2)* h(2)));
  diagConst(diagInd(end),3) = diagConst(diagInd(end),3) + (h(j)+h(j+1))/(Re*h(N+1)*h(N+1));
  %rhs
  rhsConst(diagInd(1)     ) = (h(j)+h(j+1))/(Re*h(2)*h(2)    );
  rhsConst(diagInd(end)   ) = (h(j)+h(j+1))/(Re*h(N+1)*h(N+1));
end

%strat loop to find ss solution
while (udiff>ssTol || vdiff>ssTol)
    %increment iter
    iter=iter+1;
    
    %------------------------------------------------------------------------------%
    % u component -----------------------------------------------------------------%
    %------------------------------------------------------------------------------%
    for j = 2: N
        
        diagInd = ((j-1)*(N + 2) + 2) : ((j-1)*(N + 2) + (N+1));
        
        np_interface   = 0.5 * (uold(2:N+1,j)   + uold(2:N+1,j+1));
        n_interface    = 0.5 * (uold(2:N+1,j-1) + uold(2:N+1,j));
        tp_interface   = 0.5 * (vold(2:N+1,j)   + vold(2:N+1,j+1));
        t_interface    = 0.5 * (vold(1:N,j)     + vold(1:N,j+1));
        
        if (upwind)
        np_upwind        = np_interface > 0;
        n_upwind         = n_interface >  0;
        tp_upwind        = tp_interface > 0;
        t_upwind         = t_interface >  0;
        end
        %All terms by velocity components
        %u(1:N,j)
        diag(diagInd-(N+2),1) = diagConst(diagInd-(N+2),1) -adv(j,2) * (2 * (1 - t_upwind)).*t_interface ;
        
        %u(2:N+1,j-1)
        diag(diagInd-1,2)     = diagConst(diagInd-1,2) - adv(j,1) * (2 * (1 - n_upwind)).*n_interface;
        
        %u(2:N+1,j)
        diag(diagInd,3)       = diagConst(diagInd,3) +  adv(j,2) * ...
            ((2 * np_upwind) .*np_interface ... 
             -(2 * n_upwind) .*n_interface  ...
             +(2 * tp_upwind).*tp_interface ... 
             -(2 * t_upwind) .*t_interface);
        
        %u(2:N+1,j+1)
        diag(diagInd+1,4)     = diagConst(diagInd+1,4) + adv(j,3) * (2 * (1 - np_upwind)).* np_interface;
        
        %u(3:N+2,j)
        diag(diagInd+(N+2),5) = diagConst(diagInd+(N+2),5) + adv(j,2) *(2 * (1 - tp_upwind)).* tp_interface;
        
        %RHS
        rhs(diagInd       ) = (rho*(h(j)+h(j+1))/(2 * dt)) * uold(2:N+1,j)...
            -(pres(:,j)-pres(:,j-1)).*h(2:N+1)...
            +1/Re*(vold(2:N+1,j+1)/h(j+1)-vold(2:N+1,j)/h(j))...
            -1/Re*(vold(1:N  ,j+1)/h(j+1)-vold(1:N,j  )/h(j));
        
        %Apply BC
        rhs(diagInd(1)     ) = rhs(diagInd(1)    )  + rhsConst(diagInd(1)     ) *uLO(j);
        rhs(diagInd(end)   ) = rhs(diagInd(end)  )  + rhsConst(diagInd(end)   ) *uUP(j);
    end
    
    vmat  = spdiags(diag,diagPos,(N+2)*(N+1),(N+2)*(N+1));
    if (implicit)
        %Implicit Euler
        vsol  = vmat\rhs;
        %re-map velocity u
        for j=2:N
            diagInd = ((j-1)*(N + 2) + 2) : ((j-1)*(N + 2) + (N+1));
            u(2:N+1,j) = vsol(diagInd);
        end
    else
        %Explicit Euler
        vsol(:)= initval;
        for j=2:N
            diagInd = ((j-1)*(N + 2) + 2) : ((j-1)*(N + 2) + (N+1));
            vsol(diagInd) = uold(2:N+1,j);
        end
        vexp = vmat*vsol;
        for j = 2: N
            diagInd = ((j-1)*(N + 2) + 2) : ((j-1)*(N + 2) + (N+1));
            Ru(:,j) = vexp(diagInd)...
                -1/Re*(vold(2:N+1,j+1)/h(j+1)-vold(2:N+1,j)/h(j)) +1/Re*(vold(1:N  ,j+1)/h(j+1)-vold(1:N,j  )/h(j));
            %apply BC
            Ru(1,j)    =  Ru(1,j  ) - uLO(j)*(h(j)+h(j+1))/(Re*h(2)  *h(2)  );
            Ru(end,j)  =  Ru(end,j) - uUP(j)*(h(j)+h(j+1))/(Re*h(N+1)*h(N+1));
            Ru(:,j)    =  2 * dt/(rho*(h(j)+h(j+1))) * Ru(:,j);
        end
    end
    if (debug == 1) 
        Ru
        break;
    end
    %------------------------------------------------------------------------------%
    % v component -----------------------------------------------------------------%
    %------------------------------------------------------------------------------%
    for j = 2: N
        
        diagInd = ((j-1)*(N + 2) + 2) : ((j-1)*(N + 2) + (N+1));
        
        np_interface   = 0.5 * (vold(j,2:N+1)   + vold(j+1,2:N+1))';
        n_interface    = 0.5 * (vold(j-1,2:N+1) + vold(j,2:N+1))';
        tp_interface   = 0.5 * (uold(j,2:N+1)   + uold(j+1,2:N+1))';
        t_interface    = 0.5 * (uold(j,1:N)     + uold(j+1,1:N))';
        
        if (upwind)
        np_upwind        = np_interface > 0;
        n_upwind         = n_interface >  0;
        tp_upwind        = tp_interface > 0;
        t_upwind         = t_interface >  0;
        end
        %All terms by velocity components
        %u(1:N,j)
        diag(diagInd-(N+2),1) = diagConst(diagInd-(N+2),1) -adv(j,2) * (2 * (1 - t_upwind)).*t_interface ;
        
        %u(2:N+1,j-1)
        diag(diagInd-1,2)     = diagConst(diagInd-1,2) - adv(j,1) * (2 * (1 - n_upwind)).*n_interface;
        
        %u(2:N+1,j)
        diag(diagInd,3)       = diagConst(diagInd,3) +  adv(j,2) * ...
            ((2 * np_upwind) .*np_interface ... 
             -(2 * n_upwind) .*n_interface  ...
             +(2 * tp_upwind).*tp_interface ... 
             -(2 * t_upwind) .*t_interface);
        
        %u(2:N+1,j+1)
        diag(diagInd+1,4)     = diagConst(diagInd+1,4) + adv(j,3) * (2 * (1 - np_upwind)).* np_interface;
        
        %u(3:N+2,j)
        diag(diagInd+(N+2),5) = diagConst(diagInd+(N+2),5) + adv(j,2) *(2 * (1 - tp_upwind)).* tp_interface;
        
        %RHS
        rhs(diagInd       ) = (rho*(h(j)+h(j+1))/(2 * dt)) * vold(j,2:N+1)'...
            -(pres(j,:)-pres(j-1,:))'.*h(2:N+1)...
            +1/Re*(uold(j+1,2:N+1)/h(j+1)-uold(j,2:N+1)/h(j))'...
            -1/Re*(uold(j+1,1:N  )/h(j+1)-uold(j,1:N  )/h(j))';
        
        %Apply BC
        rhs(diagInd(1)  )    = rhs(diagInd(1)    )  + rhsConst(diagInd(1)     )*vLF(j);
        rhs(diagInd(end)  )  = rhs(diagInd(end)  )  + rhsConst(diagInd(end)   )*vRH(j);
    end
    
    
    vmat  = spdiags(diag,diagPos,(N+2)*(N+1),(N+2)*(N+1));
    if (implicit)
        %Implicit Euler
        vsol  = vmat\rhs;
        %re-map velocity u
        for j=2:N
            diagInd = ((j-1)*(N + 2) + 2) : ((j-1)*(N + 2) + (N+1));
            v(j,2:N+1) = vsol(diagInd);
        end
    else
        %Explicit Euler
        vsol(:)=initval;
        for j=2:N
            diagInd = ((j-1)*(N + 2) + 2) : ((j-1)*(N + 2) + (N+1));
            vsol(diagInd) = vold(j,2:N+1);
        end
        vexp = vmat*vsol;
        for j = 2: N
            diagInd = ((j-1)*(N + 2) + 2) : ((j-1)*(N + 2) + (N+1));
            
            Rv(j,:) = vexp(diagInd)'...
                -1/Re*(uold(j+1,2:N+1)/h(j+1)-uold(j,2:N+1)/h(j))+1/Re*(uold(j+1,1:N)/h(j+1)-uold(j, 1:N )/h(j));
            %apply BC
            Rv(j,1)    =  Rv(j,1  ) - vLF(j)*(h(j)+h(j+1))/(Re*h(2)  *h(2)  );
            Rv(j,end)  =  Rv(j,end) - vRH(j)*(h(j)+h(j+1))/(Re*h(N+1)*h(N+1));
            Rv(j,:)    =   2 * dt/(rho*(h(j)+h(j+1))) * Rv(j,:);
        end
    end
    if (debug == 1) 
        Rv
        break;
    end
    
    if (~simple)
        %--------------------------------------------------------------------------%
        %-Explicit-----------------------------------------------------------------%
        %--------------------------------------------------------------------------%
        fm = (Ru(:,2:N+1) - Ru(:,1:N)) + (Rv(2:N+1,:) - Rv(1:N,:));
        f= zeros(N,1);
        for i=1:N
            for j=1:N
                k = (i-1)*N + j;
                f(k) = fm(i,j);
            end
        end
        
        % %  The folloing statment solves the Poisson equation for the pressure
        p= A\f;
        
        % %   In the next lines the presssure solution is placed in the computational domain
        
        for i=1:N
            for j=1:N
                k = (i-1)*N + j;
                pres(i,j)= p(k);
            end
        end
        
        % %
        % %  Having obtained the pressure we can update the velocity field:
        % %
        for j = 2: N
            u(2:N+1,j) =  - 2 * dt/(rho*(h(j)+h(j+1))) *(pres(:,j)-pres(:,j-1)).*h(2:N+1)  - Ru(:,j)  ;
        end
        for i = 2: N
            v(i,2:N+1) =  - 2 * dt/(rho*(h(i)+h(i+1))) *(pres(i,:)-pres(i-1,:)).*h(2:N+1)' - Rv(i,:)  ;
        end
        
    else
        %--------------------------------------------------------------------------%
        %-Simple: Semi Implicit (predictor-corrector)------------------------------%
        %--------------------------------------------------------------------------%
        diver = (u(2:N+1,2:N+1) -u(2:N+1,1:N)) + (v(2:N+1,2:N+1) - v(1:N,2:N+1));
        
        for i=1:N
            for j=1:N
                k = (i-1)*N + j;
                f(k) = -diver(i,j);
            end
        end
        
        pcvec= A\f;
        
        for i=1:N
            for j=1:N
                k = (i-1)*N + j;
                pc(i,j)= pcvec(k);
            end
        end
        
        % step 4: update the velocities with the pressure correction
        for j = 2: N
            u(2:N+1,j) =  u(2:N+1,j) -  2 * dt/(rho*(h(j)+h(j+1))) *(pc(:,j)-pc(:,j-1)).*h(2:N+1);
        end
        for i = 2: N
            v(i,2:N+1) =  v(i,2:N+1) -  2 * dt/(rho*(h(i)+h(i+1))) *(pc(i,:)-pc(i-1,:)).*h(2:N+1)';
        end
        pres = pres + relaxFac * pc;
        
    end
    
    %Check convergenge
    %check 3
    diver = (u(2:N+1,2:N+1) -u(2:N+1,1:N)) + (v(2:N+1,2:N+1) - v(1:N,2:N+1));
    maxdiv = max(max(abs(diver)));
    
    %check if we need to exit the loop
    udiff = max(max(abs(u-uold)))/(dt*(max(max(abs(u)))));
    vdiff = max(max(abs(v-vold)))/(dt*(max(max(abs(v)))));
    
    %exit if needed
    if (maxdiv>divTol)
        ['exiting bcs max diver is ' num2str(maxdiv) ' at iter ' num2str(iter) ' and rhs sum is ' num2str(sum(f)) ]
        break;
    elseif (iter>maxIter)
        [' max iterations reached!' ]
        break;
    end
    %update for next iteration
    uold=u;
    vold=v;
end

%print iteration
['iter ' num2str(iter) ' divergence ' num2str(maxdiv)]

%plot velocity vectors
 xc = zeros(N,N);
 yc = zeros(N,N);
 uc = zeros(N,N);
 vc = zeros(N,N);
 for i=1:N
     xc(:,i) = x(i)+h(i+1)/2;
     yc(i,:) = x(i)+h(i+1)/2;
     uc(:,i) = (u(2:N+1,i) + u(2:N+1,i+1))./2;
     vc(i,:) = (v(i,2:N+1) + v(i+1,2:N+1))./2;
 end
 quiver(xc,yc,uc,vc);

%plot reference (ghia et al. 1982, High-Re Solutions for Incompressible Flow
%Using the Navier-Stokes Equations and a
%Multigrid Method, pg 398): velocity values (+) and obtained solution
figure
subplot(1,2,1); 
plot(u(2:N+1,N/2+1),x(1:end-1));
hold on
ugRe100=[ 
0.06250 -0.04192, 
0.1016  -0.06434, 
0.2813  -0.15662,
0.5000  -0.20581,
0.7344   0.00332,
0.9531   0.68717, 
0.8516   0.2315, 
0.9688   0.78871];
plot(ugRe100(:,2),ugRe100(:,1),'+r');
title('u (+ reference N=128)');xlabel('u'); ylabel('y');


subplot(1,2,2); 
plot(v(N/2+1,2:N+1),x(1:end-1));
hold on
ygRe100=[ 
0.06250   0.09233,
0.0781    0.10890,
0.1563    0.16077,
0.2344    0.17527,
0.5000    0.05454, 
0.8047   -0.24533,
0.9063   -0.16914,
0.9531   -0.08864,
0.9688   -0.05906
];
plot(ygRe100(:,2),ygRe100(:,1),'+r');

title('v (+ reference N=128)');xlabel('v'); ylabel('x');