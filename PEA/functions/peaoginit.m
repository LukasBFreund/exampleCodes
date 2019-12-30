function theta=peaoginit(e,param);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %
%                 Parameters of the algorithm                 %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncont	= 3;              		% # of static equations
nbend	= 1;               		% # endogenous predetermined variables
nshoc	= 1;               		% # of shocks
nback	= nbend+nshoc;				% # of state variables
nforw	= 1;               		% # of costate variables
nstat	= nback+nforw;     		% # of state and costate variables
Mcc	= zeros(ncont,ncont);
Mcs	= zeros(ncont,nstat);
Mss0	= zeros(nstat,nstat);
Mss1	= zeros(nstat,nstat);
Msc0	= zeros(nstat,ncont);
Msc1	= zeros(nstat,ncont);
Mse	= zeros(nstat,nshoc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %
%             structural parameters of the economy            %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ab		= param(1); 
alpha = param(2);
beta 	= param(3);
delta	= param(4);
rho	= param(5);
se		= param(6);
sigma	= param(7);
long	= param(8);
init	= param(9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %
%                          steady state                       %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ysk =(1-beta*(1-delta))/(alpha*beta);
ksy = 1/ysk;
ys	 = ksy^(alpha/(1-alpha));
ks  = ys^(1/alpha);
is  = delta*ks;
cs  = ys-is;
ls	 = cs^(-sigma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %
%                      matrix coefficients                    %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% y i 
% 1 2 
%
% k z c
% 1 2 3
%
% Output
%
Mcc(1,1)	= 1;
Mcs(1,1)	= alpha;
Mcs(1,2)	= 1;
%
% investment
%
Mcc(2,1)	= 1;
Mcc(2,2)	= -is/ys;
Mcc(2,3)	= -cs/ys;
%
% consumption
%
Mcc(3,3)	= -sigma;
Mcs(3,3)	= 1;
%
% capital
%
Mss0(1,1)	=	1;
Mss1(1,1)	=	delta-1;
Msc1(1,2)	=	delta;
%
% technology shock
%
Mss0(2,2)	= 1;
Mss1(2,2)	= -rho;
Mse(2,1)		= 1;
%
% Euler
%
Mss0(3,1)	= (1-beta*(1-delta));
Mss0(3,3)	= -1;
Mss1(3,3)	= 1;
Msc0(3,1)	= (1-beta*(1-delta));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %
%            Solving the system                               %
%                                                             %
%           -1                                                %
% X(t) = Mcc  Mcs S(t)                                        %
%                                                             %
%                      -1             -1                      %
% S(t+1)=(Mss0-Msc0 Mcc  Mcs)(Msc1 Mcc  Mcs-Mss1)S(t)         %
%                       -1                                    %
%        +(Mss0-Msc0 Mcc  Mcs)e(t+1)                          %
%                                                             %
% S(t+1) = W S(t) + R e(t+1)                                  %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M0=inv(Mss0-Msc0*inv(Mcc)*Mcs);
M1=(Mss1-Msc1*inv(Mcc)*Mcs);
W=-M0*M1;
%
% MU -> eigenvalues, P -> eigenvectors
%
[P,MU] = eig(W);
AMU=abs(MU);
[MU,k] = sort(diag(AMU));
P=P(:,k);
Q=inv(P);
%
% Direct solution
%
Gamma=-inv(Q(nback+1:nstat,nback+1:nstat))*Q(nback+1:nstat,1:nback);
MSS=W(1:nback,1:nback)+W(1:nback,nback+1:nstat)*Gamma;
PI=inv(Mcc)*(Mcs(:,1:nback)+Mcs(:,nback+1:nstat)*Gamma);
MSE=[zeros(nbend,nshoc);eye(nshoc)];

S			= zeros(nback,long+init);
S(:,1)	= MSE*e(1);
for i		= 2:long+init;
   S(:,i)= MSS*S(:,i-1)+MSE*e(i);
end;
lb			= Gamma*S;
lb			= ls*exp(lb);lb=lb(:);
k			= log(ks)+S(1,:);k=k(:);
ek			= exp(k);
a			= S(2,:);a=a(:);
ea			= exp(a);
T			= init+1:init+long-1;
T1			= init+2:init+long;
X			= [ones(long-1,1) k(T) a(T) k(T).*k(T) a(T).*a(T) k(T).*a(T)];
y			= log(beta*lb(T1).*(alpha*ea(T1).*ek(T1).^(alpha-1)+1-delta));
theta		= X\y;
