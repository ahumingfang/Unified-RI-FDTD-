clc;clear;close all

NZ=3000; %computational domain
Nt=10000; %time steps

% matrices
EZ  =zeros(2,NZ);
HY  =zeros(1,NZ);

GEZ =zeros(1,NZ);
GP =zeros(1,NZ);
Hyx=zeros(1,NZ);

m0=12.566e-7; %permeability of fREe_xe space
e0=8.854187817620e-12; %permittivity of fREe_xe space
c0=3.0e8;

dx=0.1e-3;
dt=1*dx/(c0); %time step

ki=0.1e-9;
for n=1:Nt
    func(n)=100*exp(-( (n*dt-5.0*ki)/(ki) )^2);
end

% thichness
za=NZ/2;
zb=za+1;

%the source location
zsource=za-1000;

% parameters
e0_inf=29.9;
es=47.9;
tao=43.6e-12;
deltae=es-e0_inf;
sigma_s=0.540;

K=e0_inf;
sigma=deltae;
Alfa=1;
tao1=tao;

CA=zeros(1,NZ-1);
CA(1:za-1)=1;
CA(za:zb)=(e0/dt-sigma_s/2)/(e0/dt+sigma_s/2);
CA(zb+1:NZ-1)=1;

CB=zeros(1,NZ-1);
CB(1:za-1)=dt/e0;
CB(za:zb)=1/(e0/dt+sigma_s/2);
CB(zb+1:NZ-1)=dt/e0;

CQ=dt/m0;
% PML Parameters
alpha_1 = 4.0; 
alpha_1_aa = 1.0; 
 
% PML thickness
NPML = 11; 

max_sigma_x = 1.0 * (alpha_1 + 1) / (150 * pi * dx);
max_alpha_1_x = 0.005; 
%% Initialization
sigma_x = zeros(1,NZ);
alpha_1_x = zeros(1,NZ);
sigma_mx = zeros(1,NZ);
alpha_1_mx = zeros(1,NZ);

RAe_x = zeros(1,NZ);
RBe_x = zeros(1,NZ);
RCe_x = zeros(1,NZ);
RDe_x = zeros(1,NZ);
REe_x = zeros(1,NZ);
RFe_x = zeros(1,NZ);
sigma_s_index = zeros(1,NZ);

RAe_x(za:zb)=((2*tao1+dt*Alfa)/(2*tao1*K+dt*(Alfa*K+sigma)))-1;
RBe_x(za:zb)=(2*tao1*K)/(2*tao1*K+dt*(Alfa*K+sigma));
RCe_x(za:zb)=dt*(Alfa-Alfa*K-sigma)/(tao1*K);
RDe_x(za:zb)=dt*(Alfa*K+sigma)/(tao1*K);
REe_x(za:zb)=1-RDe_x(za:zb).*RBe_x(za:zb);
RFe_x(za:zb)=RCe_x(za:zb)-RAe_x(za:zb).*RDe_x(za:zb);
sigma_s_index(za:zb)=sigma_s;
%% E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at front boundary (x)
for i = 1:NPML
    sigma_x(i) = max_sigma_x * (((NPML) - i) / (NPML - 1.0))^alpha_1;
	alpha_1_x(i) = max_alpha_1_x * (((NPML) - i) / (NPML - 1.0))^alpha_1_aa;
    RAe_x(i) = (2.0 * e0 + dt * alpha_1_x(i)) / (2.0 * e0 + dt * (alpha_1_x(i) + sigma_x(i))) - 1.0;
	RBe_x(i) = (2.0 * e0) / (2.0 * e0 + dt * (alpha_1_x(i) + sigma_x(i)));
	RCe_x(i) = dt * (alpha_1_x(i) - alpha_1_x(i) * 1.0 - sigma_x(i)) / e0;
	RDe_x(i) = dt * (alpha_1_x(i) + sigma_x(i)) / e0;
	REe_x(i) = 1.0 - RDe_x(i) * RBe_x(i);
	RFe_x(i) = RCe_x(i) - RAe_x(i) * RDe_x(i);
end

% at back boundary (x)
for i = NZ - NPML + 1: NZ
	sigma_x(i) = sigma_x(NZ + 1 - i);
	alpha_1_x(i) = alpha_1_x(NZ + 1 - i);
	RAe_x(i) = (2.0 * e0 + dt * alpha_1_x(i)) / (2.0 * e0 + dt * (alpha_1_x(i) + sigma_x(i))) - 1.0;
	RBe_x(i) = (2.0 * e0) / (2.0 * e0 + dt * (alpha_1_x(i) + sigma_x(i)));
	RCe_x(i) = dt * (alpha_1_x(i) - alpha_1_x(i) * 1.0 - sigma_x(i)) / e0;
	RDe_x(i) = dt * (alpha_1_x(i) + sigma_x(i)) / e0;
	REe_x(i) = 1.0 - RDe_x(i) * RBe_x(i);
	RFe_x(i) = RCe_x(i) - RAe_x(i) * RDe_x(i);
end

RAh_x = zeros(1,NZ);
RBh_x = zeros(1,NZ);
RCh_x = zeros(1,NZ);
RDh_x = zeros(1,NZ);
REh_x = zeros(1,NZ);
RFh_x = zeros(1,NZ);
%%  H  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at front boundary (x)
for i = 1:NPML-1
    sigma_mx(i) = max_sigma_x * (((NPML) - i - 0.5) / (NPML - 1.0))^alpha_1;
	alpha_1_mx(i) = max_alpha_1_x * (((NPML) - i - 0.5) / (NPML - 1.0))^alpha_1_aa;
	RAh_x(i) = (2.0 * e0 + dt * alpha_1_mx(i)) / (2.0 * e0 + dt * (alpha_1_mx(i) + sigma_mx(i))) - 1.0;
	RBh_x(i) = (2.0 * e0) / (2.0 * e0 + dt * (alpha_1_mx(i) + sigma_mx(i)));
	RCh_x(i) = dt * (alpha_1_mx(i) - alpha_1_mx(i) * 1.0 - sigma_mx(i)) / e0;
	RDh_x(i) = dt * (alpha_1_mx(i) + sigma_mx(i)) / e0;
	REh_x(i) = 1.0 - RDh_x(i) * RBh_x(i);
	RFh_x(i) = RCh_x(i) - RAh_x(i) * RDh_x(i);
end

% at back boundary (x)
for i = NZ - NPML + 1 :  NZ - 1
	sigma_mx(i) = sigma_mx(NZ - i);
	alpha_1_mx(i) = alpha_1_mx(NZ - i);
	RAh_x(i) = (2.0 * e0 + dt * alpha_1_mx(i)) / (2.0 * e0 + dt * (alpha_1_mx(i) + sigma_mx(i))) - 1.0;
	RBh_x(i) = (2.0 * e0) / (2.0 * e0 + dt * (alpha_1_mx(i) + sigma_mx(i)));
	RCh_x(i) = dt * (alpha_1_mx(i) - alpha_1_mx(i) * 1.0 - sigma_mx(i)) / e0;
	RDh_x(i) = dt * (alpha_1_mx(i) + sigma_mx(i)) / e0;
	REh_x(i) = 1.0 - RDh_x(i) * RBh_x(i);
	RFh_x(i) = RCh_x(i) - RAh_x(i) * RDh_x(i);
end
%%
sigma_s_index_CB_RA = zeros(1,NZ);
sigma_s_index_CB_RA(2:NZ-1)=sigma_s_index(2:NZ-1) .* CB(2:NZ-1) .* RAe_x(2:NZ-1);
sigma_s_index_CB_RB = zeros(1,NZ);
sigma_s_index_CB_RB(2:NZ-1)=sigma_s_index(2:NZ-1) .* CB(2:NZ-1) .* RBe_x(2:NZ-1);

RI_EZ_trans_1=zeros(1,Nt);

for n=1:Nt
    n

    EZ(2,:)= EZ(1,:); 
    EZ(1,2:NZ-1) = CA(2:NZ-1).*EZ(2,2:NZ-1) + CB(2:NZ-1).*(HY(1,2:NZ-1)-HY(1,1:NZ-2))/dx + ...
                   CB(2:NZ-1).*RAe_x(2:NZ-1).*(HY(1,2:NZ-1)-HY(1,1:NZ-2))/dx + CB(2:NZ-1).* RBe_x(2:NZ-1).* GEZ(2:NZ-1) - ...
                   sigma_s_index_CB_RA(2:NZ-1).* EZ(2,2:NZ-1) -  sigma_s_index_CB_RB(2:NZ-1).* GP(2:NZ-1);
      
   % Auxiliary variable          
   GP(2:NZ-1) =  REe_x(2:NZ-1) .* GP(2:NZ-1)  + RFe_x(2:NZ-1) .* EZ(2,2:NZ-1);
   GEZ(2:NZ-1) = REe_x(2:NZ-1) .* GEZ(2:NZ-1) + RFe_x(2:NZ-1) .* ((HY(1,2:NZ-1)-HY(1,1:NZ-2))/dx); 

   EZ(1,zsource)=EZ(1,zsource)+func(n);
   RI_EZ_trans_1(n)=EZ(1,zb+12);
   
   % H component
   HY(1,1:NZ-1)=HY(1,1:NZ-1) + CQ * (EZ(1,2:NZ)-EZ(1,1:NZ-1))/dx+ CQ * (RAh_x(1:NZ-1) .*(EZ(1,2:NZ)-EZ(1,1:NZ-1))/dx + RBh_x(1:NZ-1) .* Hyx(1:NZ-1) );
   
   % Auxiliary variable 
   Hyx(1:NZ-1) = REh_x(1:NZ-1) .* Hyx(1:NZ-1) + RFh_x(1:NZ-1) .* (EZ(1,2:NZ)-EZ(1,1:NZ-1))/dx;         
end

