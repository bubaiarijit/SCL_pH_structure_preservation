clc
clear all
close all

%%
load('mechanical_example_data.mat');
clear n
n =40;
save('mechanical_example_data.mat')
global n
tspan=0:.01:100;
syms s t

%%

K1=K(1:n,1:n);

M1=M(1:n,1:n);
M_inv=inv(M1);
alpha=10^-20;
epsilon_o=1000   ;
epsilon_c=5000;
global beta
beta = [];
for i=1:n
    beta = [beta 0.2*10^(-10)*M1(i,i)*9.81]; 
end
G=[0;1;zeros(n-2,1)];
global M1
global K1
global G

B = [zeros(n,1); G];
H = [K1 zeros(n,n);zeros(n,n) M_inv];


%% Generalized controllability Gramian
P=sdpvar(2*n);
A=[zeros(n,n) eye(n,n);-eye(n,n) -0.1*eye(n,n)]*H;
JR=[zeros(n,n) eye(n,n);-eye(n,n) -0.1*eye(n,n)];
F=[];
F=[F, 
    A*P+P*A'+B*B'<=-10^(-3)*P
    ];
for i=1:n
    r=0.1*ones(n,1);
    r(i)= r(i)+(beta(i)/(sqrt(alpha))* M1(i,i));
    A=[zeros(n,n) eye(n,n);-eye(n,n) -diag(r)]*H;
    F=[F,
        A*P+P*A'+B*B'<=-10^(-3)*P
        ];
end
F=[F, 
    P>=10^(-10)*eye(2*n)];
options=sdpsettings('solver','sedumi');
sol=optimize(F,trace(P),options);
if sol.problem == 0
disp('Feasible')
elseif sol.problem == 1
disp('Infeasible')
else
disp('Something else happened')    
end

P=value(P);

%% Generalized observability Gramian
% Q=sdpvar(2*n);
% A=[zeros(n,n) eye(n,n);-eye(n,n) -0.1*eye(n)]*H;
% JR=[zeros(n,n) eye(n,n);-eye(n,n) -0.1*eye(n)];
% 
% F=[];
% F=[F,
%     Q*A+A'*Q+H*B*B'*H<=-10^(-3)*Q
%     ];
% for i=1:n
%     r=0.1*ones(n,1);
%     r(i)= r(i)+(beta(i)/(sqrt(alpha))* M1(i,i));
%     A=[zeros(n,n) eye(n,n);-eye(n,n) -diag(r)]*H;
%     F=[F,
%          Q*A+A'*Q+H*B*B'*H<=-10^(-3)*Q
%          ];
% end
% F=[F, 
%     Q>=10^(-10)*eye(2*n)];
% options=sdpsettings('solver','sedumi');
% sol=optimize(F,trace(Q),options);
% if sol.problem == 0
% disp('Feasible')
% elseif sol.problem == 1
% disp('Infeasible')
% else
% disp('Something else happened')    
% end
% 
% Q=value(Q);

%% ---------Structure-preserving generalized balancing(Algorithm 1)-------- 
 PhP=chol(P);
[UHP,LHP]=svd(PhP*H*PhP');
L2PQ=diag(sdpvar(2*n,1));
JR=[zeros(n,n) eye(n,n);-eye(n,n) -0.1*eye(n,n)];
JRc=UHP'*inv(PhP')*JR*inv(PhP)*UHP;
Bc=UHP'*inv(PhP')*B
F=[];
F=[F, 
    L2PQ*JRc*LHP+LHP*JRc'*L2PQ+Bc*Bc'<=-10^(-3)*eye(2*n)
    ];
for i=1:n
    r=0.1*ones(n,1);
    r(i)= r(i)+(beta(i)/(sqrt(alpha))* M1(i,i));
    JR=[zeros(n,n) eye(n,n);-eye(n,n) -diag(r)];
    JRc=UHP'*inv(PhP')*JR*inv(PhP)*UHP;
    F=[F,
        L2PQ*JRc*LHP+LHP*JRc'*L2PQ+Bc*Bc'<=-10^(-3)*eye(2*n)
        ];
end
F=[F, 
    L2PQ>=10^(-10)*eye(2*n)];
sol=optimize(F,trace(L2PQ));
if sol.problem == 0
disp('Feasible')
elseif sol.problem == 1
disp('Infeasible')
else
disp('Something else happened')    
end

L2PQ=value(L2PQ);
L2PQ=diag(sort(diag(L2PQ),'descend'));
LPQ=sqrt(L2PQ);
Qsp=inv(PhP)*UHP*L2PQ*UHP'*inv(PhP');
global Wsp
Wsp=PhP'*UHP*sqrt(inv(LPQ));
Hgc=Wsp'*H*Wsp;
Qgc=Wsp'*Qsp*Wsp;
Pgc=inv(Wsp)*P*inv(Wsp');
%% ----Structure-preserving generalized balancing(Algorithm 2)----
% PhQ=chol(Q);
% [UHQ,LHQ]=svd(inv(PhQ')*H*inv(PhQ));
% L2PQ=diag(sdpvar(2*n,1));
% JR=[zeros(n,n) eye(n,n);-eye(n,n) -0.1*eye(n,n)];
% JRo=UHQ'*PhQ*JR*PhQ'*UHQ;
% Bo=UHQ'*PhQ*B;
% F=[];
% F=[F, 
%    JRo*LHQ*L2PQ+L2PQ*LHQ*JRo'+Bo*Bo'<=-10^(-3)*eye(2*n)
%     ];
% for i=1:n
%     r=0.1*ones(n,1);
%     r(i)= r(i)+(beta(i)/(sqrt(alpha))* M1(i,i));
%     JR=[zeros(n,n) eye(n,n);-eye(n,n) -diag(r)];
%     JRo=UHQ'*PhQ*JR*PhQ'*UHQ;
%     F=[F,
%         JRo*LHQ*L2PQ+L2PQ*LHQ*JRo'+Bo*Bo'<=-10^(-3)*eye(2*n)
%         ];
% end
% F=[F, 
%     L2PQ>=10^(-10)*eye(2*n)];
% sol=optimize(F,trace(L2PQ));
% if sol.problem == 0
% disp('Feasible')
% elseif sol.problem == 1
% disp('Infeasible')
% else
% disp('Something else happened')    
% end
% 
% L2PQ=value(L2PQ);
% L2PQ=diag(sort(diag(L2PQ),'descend'));
% LPQ=sqrt(L2PQ);
% Psp=inv(PhQ)*UHQ*L2PQ*UHQ'*inv(PhQ');
% global Wsp
% Wsp=inv(PhQ)*UHQ*sqrt(LPQ);
% Hgc=Wsp'*H*Wsp;
% Qgc=Wsp'*Q*Wsp;
% Pgc=inv(Wsp)*Psp*inv(Wsp');

%% --------- Figure 1 ---------
z0=zeros(2*n,1);
Nr1 = 30;
Nr2 = 20;
Cbal2=B'*H*Wsp;
[t,zeg]=ode45(@(t,z) func(t,z,Wsp), tspan, z0);
yeg=[];
for i=1:length(zeg(:,1))
    yeg=[yeg;Cbal2*zeg(i,:)'];
end
figure(1)
plot(t,yeg,'Linewidth',3)
set(gca,'FontSize',40)
grid on
xlabel('time', 'fontsize',50)
ylabel('output', 'fontsize',50)
hold on
legend('Original System','FontSize',30)
hold on
[t,zreg]=ode45(@(t,z) funcred(t,z,Wsp,Nr1), tspan, z0);
yreg=[];
for i=1:length(zreg(:,1))
    yreg=[yreg;Cbal2*zreg(i,:)'];
end    
plot(t,yreg,'--','DisplayName',strcat('Reduced system, dim=',num2str(Nr1)),'Linewidth',3)
hold on
[t,zreg]=ode45(@(t,z) funcred(t,z,Wsp,Nr2), tspan, z0);
yreg=[];
for i=1:length(zreg(:,1))
    yreg=[yreg;Cbal2*zreg(i,:)'];
end    
plot(t,yreg,'-.','DisplayName',strcat('Reduced system, dim=',num2str(Nr2)),'Linewidth',3)
hold off
%% Figure 2
figure(2)
d = diag(LPQ);
plot(d,'--*','Linewidth',3)
set(gca,'FontSize',40)
grid on
title('The diagonal entries of \Lambda_{PQ}')
xlabel('State number (i)')
ylabel('\sigma_i value')

%% Functions 
function zdot = func(t,z,W)
alpha=10^-20;
global n
global K1
global M1
global beta
M_inv=inv(M1);
%global R
global G
    u = sin(t);
    %u=10*exp(-5*t);
    q = z(1:n);
    p = z(n:2*n);
    x = W*z;
    
    r=[];
    for k = n+1:2*n
        r = [r beta(k-n)/sqrt(alpha+x(k)^2/M_inv(k-n,k-n)^2)];
    end
%     for k = n+1:2*n
%         r = [r sign(z(k))/M];
%     end    
    R = diag(r);
    zdot = inv(W)*[zeros(n,n) eye(n,n);-eye(n,n) -(0.1*eye(n)+R)]*[K1 zeros(n,n);zeros(n,n) M_inv]*W*z + inv(W)*[zeros(n,1); G]*u;
end  

function zrdot = funcred(t,z,W,Nr)
alpha=10^-20;
global n
global K1
global M1
global beta
M_inv=inv(M1);
%global R
global G
    u = sin(t);
    %u=10*exp(-5*t);
    z(Nr+1:2*n)=zeros(2*n-Nr,1);
    x = W*z;
    
    r=[];
    for k = n+1:2*n
        r = [r beta(k-n)/sqrt(alpha+x(k)^2/M_inv(k-n,k-n)^2)];
    end
%     for k = n+1:2*n
%         r = [r sign(z(k))/M];
%     end    
    R = diag(r);
    Abal=inv(W)*[zeros(n,n) eye(n,n);-eye(n,n) -(0.1*eye(n)+R)]*[K1 zeros(n,n);zeros(n,n) M_inv]*W;
    Abal(Nr+1:2*n,:)=zeros(2*n-Nr,2*n);
    Ar=Abal;
    Bbal=inv(W)*[zeros(n,1); G];
    Bbal(Nr+1:2*n)=zeros(2*n-Nr,1);
    Br=Bbal;
    zrdot = Ar*z + Br*u;
    %zdot = [zeros(n,n) eye(n,n);-eye(n,n) zeros(n,n)]*[K zeros(n,n);zeros(n,n) M_inv]*z - [zeros(n,1); r']+ [zeros(n,1); G]*u;
end  
function zdot = funclinear(t,z)
global n
global K1
global M1
M_inv=inv(M1);
%global R
global G
    u = 20*sin(t);
    zdot = [zeros(n,n) eye(n,n);-eye(n,n) -0.1*eye(n)]*[K1 zeros(n,n);zeros(n,n) M_inv]*z + [zeros(n,1); G]*u;
end  