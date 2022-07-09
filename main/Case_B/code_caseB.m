clc;
clear; 
close all;  
 
% Grid for EM fields
Imax = 191;  
Jmax = 351; 
Kmax = 191; 

Ntimesteps =3824; 
N=20;    
epsR = 1.0; 
sigM1 = 0.0;  
 
pi =  3.1415926535897932384626433832795; 
c0 = 2.99792458E8; 
mu0 = 4.0 * pi * 1.0E-7;  
eps0 = 1.0/(c0*c0*mu0); 

MUR0 = 4.0 * pi * 1.0E-7;            % permeability
EPSILON0 = 1.0 / (c0 * c0 * mu0);    % permittivity

% Cell Size
dx = 0.2E-3;
dy = 0.2E-3; 
dz = 0.2E-3; 

% CFL Limit
CFL=1.0;
%dt = CD*dx /(3^(1/2)*3e8); 
% EM_DELTAT = 1 * CFL * 1 / c0 / (((1 /dx/dx + 1 /dy/dy + 1 /dz/dz))^0.5); 
%        dt = 1 * CFL * 1 / c0 / (((1 /dx/dx + 1 /dy/dy + 1 /dz/dz))^0.5); 

EM_DELTAT = 0.34e-12;
       dt = 0.34e-12;

SPREAD = 1.0E-10;                 % pulse width  ***  5.0*1.e-15*ref_s (narrow)  15.e-16
T0 = 4.0 * SPREAD;                % pulse center ***

% Computational Center            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isource01 = round(Imax/2);  
jsource01 = round(Jmax/2); 
ksource01 = round(Kmax/2); 

% PML thickness in each direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NPML = 21; 

% PML Parameters
alpha = 4.0; 
alpha_aa = 1.0; 
 
% max_sigma_x = (0.8*(alpha+1)/(dx*(muO/epsO*epsR)^0.5)); 
% max_sigma_y = (0.8*(alpha+1)/(dy*(muO/epsO*epsR)^0.5)); 
% max_sigma_z = (0.8*(alpha+1)/(dz*(muO/epsO*epsR)^0.5)); 

max_sigma_x = 1.0 * (alpha + 1) / (150 * pi * dx);
max_sigma_y = 1.0 * (alpha + 1) / (150 * pi * dy);
max_sigma_z = 1.0 * (alpha + 1) / (150 * pi * dz); 

max_alpha_x = 0.005; 
max_alpha_y = 0.005; 
max_alpha_z = 0.005;

% kappa_x_max = 1; 
% kappa_y_max = kappa_x_max;  
% kappa_z_max = kappa_x_max; 
 
% Field Components
Ex =zeros(Imax+1, Jmax+1, Kmax+1);
Ey =zeros(Imax+1, Jmax+1, Kmax+1);
Ez =zeros(Imax+1, Jmax+1, Kmax+1);
Hx =zeros(Imax+1, Jmax+1, Kmax+1);
Hy =zeros(Imax+1, Jmax+1, Kmax+1);
Hz =zeros(Imax+1, Jmax+1, Kmax+1);

Ex_1 =zeros(Imax+1, Jmax+1, Kmax+1);
Ey_1 =zeros(Imax+1, Jmax+1, Kmax+1);
Ez_1 =zeros(Imax+1, Jmax+1, Kmax+1);

Hxy=zeros(Imax+1, Jmax+1, Kmax+1);
Hxz=zeros(Imax+1, Jmax+1, Kmax+1);
Hyx=zeros(Imax+1, Jmax+1, Kmax+1);
Hyz=zeros(Imax+1, Jmax+1, Kmax+1);
Hzx=zeros(Imax+1, Jmax+1, Kmax+1);
Hzy=zeros(Imax+1, Jmax+1, Kmax+1);

CA_Tx=zeros(Imax+1, Jmax+1, Kmax+1);
CB_Tx=zeros(Imax+1, Jmax+1, Kmax+1);
CA_Ty=zeros(Imax+1, Jmax+1, Kmax+1);
CB_Ty=zeros(Imax+1, Jmax+1, Kmax+1);
CA_Tz=zeros(Imax+1, Jmax+1, Kmax+1);
CB_Tz=zeros(Imax+1, Jmax+1, Kmax+1);

GP_Sx=zeros(Imax+1, Jmax+1, Kmax+1);
GEZ_Sxy=zeros(Imax+1, Jmax+1, Kmax+1);
GEZ_Sxz=zeros(Imax+1, Jmax+1, Kmax+1);

GP_Sy=zeros(Imax+1, Jmax+1, Kmax+1);
GEZ_Syx=zeros(Imax+1, Jmax+1, Kmax+1);
GEZ_Syz=zeros(Imax+1, Jmax+1, Kmax+1);

GP_Sz=zeros(Imax+1, Jmax+1, Kmax+1);
GEZ_Szx=zeros(Imax+1, Jmax+1, Kmax+1);
GEZ_Szy=zeros(Imax+1, Jmax+1, Kmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%

CB_index_X=zeros(Imax+1, Jmax+1, Kmax+1);
RA_index_X=zeros(Imax+1, Jmax+1, Kmax+1);
RB_index_X=zeros(Imax+1, Jmax+1, Kmax+1);

RE_index_X=zeros(Imax+1, Jmax+1, Kmax+1);
RF_index_X=zeros(Imax+1, Jmax+1, Kmax+1);
sigma_s_index_X=zeros(Imax+1, Jmax+1, Kmax+1);

RA_index_X_y=zeros(Imax+1, Jmax+1, Kmax+1);
RB_index_X_y=zeros(Imax+1, Jmax+1, Kmax+1);
RE_index_X_y=zeros(Imax+1, Jmax+1, Kmax+1);
RF_index_X_y=zeros(Imax+1, Jmax+1, Kmax+1);

RA_index_X_z=zeros(Imax+1, Jmax+1, Kmax+1);
RB_index_X_z=zeros(Imax+1, Jmax+1, Kmax+1);
RE_index_X_z=zeros(Imax+1, Jmax+1, Kmax+1);
RF_index_X_z=zeros(Imax+1, Jmax+1, Kmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CB_index_Y=zeros(Imax+1, Jmax+1, Kmax+1);
RA_index_Y=zeros(Imax+1, Jmax+1, Kmax+1);
RB_index_Y=zeros(Imax+1, Jmax+1, Kmax+1);

RE_index_Y=zeros(Imax+1, Jmax+1, Kmax+1);
RF_index_Y=zeros(Imax+1, Jmax+1, Kmax+1);
sigma_s_index_Y=zeros(Imax+1, Jmax+1, Kmax+1);

RA_index_Y_x=zeros(Imax+1, Jmax+1, Kmax+1);
RB_index_Y_x=zeros(Imax+1, Jmax+1, Kmax+1);
RE_index_Y_x=zeros(Imax+1, Jmax+1, Kmax+1);
RF_index_Y_x=zeros(Imax+1, Jmax+1, Kmax+1);

RA_index_Y_z=zeros(Imax+1, Jmax+1, Kmax+1);
RB_index_Y_z=zeros(Imax+1, Jmax+1, Kmax+1);
RE_index_Y_z=zeros(Imax+1, Jmax+1, Kmax+1);
RF_index_Y_z=zeros(Imax+1, Jmax+1, Kmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CB_index_Z=zeros(Imax+1, Jmax+1, Kmax+1);
RA_index_Z=zeros(Imax+1, Jmax+1, Kmax+1);
RB_index_Z=zeros(Imax+1, Jmax+1, Kmax+1);

RE_index_Z=zeros(Imax+1, Jmax+1, Kmax+1);
RF_index_Z=zeros(Imax+1, Jmax+1, Kmax+1);
sigma_s_index_Z=zeros(Imax+1, Jmax+1, Kmax+1);

RA_index_Z_x=zeros(Imax+1, Jmax+1, Kmax+1);
RB_index_Z_x=zeros(Imax+1, Jmax+1, Kmax+1);
RE_index_Z_x=zeros(Imax+1, Jmax+1, Kmax+1);
RF_index_Z_x=zeros(Imax+1, Jmax+1, Kmax+1);

RA_index_Z_y=zeros(Imax+1, Jmax+1, Kmax+1);
RB_index_Z_y=zeros(Imax+1, Jmax+1, Kmax+1);
RE_index_Z_y=zeros(Imax+1, Jmax+1, Kmax+1);
RF_index_Z_y=zeros(Imax+1, Jmax+1, Kmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RAe_x = zeros(1,Imax+1);
RBe_x = zeros(1,Imax+1);
RCe_x = zeros(1,Imax+1);
RDe_x = zeros(1,Imax+1);
REe_x = zeros(1,Imax+1);
RFe_x = zeros(1,Imax+1);

RAe_y = zeros(1,Jmax+1);
RBe_y = zeros(1,Jmax+1);
RCe_y = zeros(1,Jmax+1);
RDe_y = zeros(1,Jmax+1);
REe_y = zeros(1,Jmax+1);
RFe_y = zeros(1,Jmax+1);

RAe_z = zeros(1,Kmax+1);
RBe_z = zeros(1,Kmax+1);
RCe_z = zeros(1,Kmax+1);
RDe_z = zeros(1,Kmax+1);
REe_z = zeros(1,Kmax+1);
RFe_z = zeros(1,Kmax+1);

RAh_x = zeros(1,Imax+1);
RBh_x = zeros(1,Imax+1);
RCh_x = zeros(1,Imax+1);
RDh_x = zeros(1,Imax+1);
REh_x = zeros(1,Imax+1);
RFh_x = zeros(1,Imax+1);

RAh_y = zeros(1,Jmax+1);
RBh_y = zeros(1,Jmax+1);
RCh_y = zeros(1,Jmax+1);
RDh_y = zeros(1,Jmax+1);
REh_y = zeros(1,Jmax+1);
RFh_y = zeros(1,Jmax+1);

RAh_z = zeros(1,Kmax+1);
RBh_z = zeros(1,Kmax+1);
RCh_z = zeros(1,Kmax+1);
RDh_z = zeros(1,Kmax+1);
REh_z = zeros(1,Kmax+1);
RFh_z = zeros(1,Kmax+1);

sigma_x = zeros(1,Imax+1);
sigma_y = zeros(1,Jmax+1);
sigma_z = zeros(1,Kmax+1);

sigma_mx = zeros(1,Imax+1);
sigma_my = zeros(1,Jmax+1);
sigma_mz = zeros(1,Kmax+1);

alpha_x = zeros(1,Imax+1);
alpha_y = zeros(1,Jmax+1);
alpha_z = zeros(1,Kmax+1);

alpha_mx = zeros(1,Imax+1);
alpha_my = zeros(1,Jmax+1);
alpha_mz = zeros(1,Kmax+1);

%% E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at front boundary (x)
for i = 1:NPML
    sigma_x(i) = max_sigma_x * (((NPML) - i) / (NPML - 1.0))^alpha;
	alpha_x(i) = max_alpha_x * (((NPML) - i) / (NPML - 1.0))^alpha_aa;
	% alpha_x(i) = max_alpha_x;
    RAe_x(i) = (2.0 * EPSILON0 + EM_DELTAT * alpha_x(i)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_x(i) + sigma_x(i))) - 1.0;
	RBe_x(i) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_x(i) + sigma_x(i)));
	RCe_x(i) = EM_DELTAT * (alpha_x(i) - alpha_x(i) * 1.0 - sigma_x(i)) / EPSILON0;
	RDe_x(i) = EM_DELTAT * (alpha_x(i) + sigma_x(i)) / EPSILON0;
	REe_x(i) = 1.0 - RDe_x(i) * RBe_x(i);
	RFe_x(i) = RCe_x(i) - RAe_x(i) * RDe_x(i);
end

% at back boundary (x)
for i = Imax - NPML + 1: Imax
	sigma_x(i) = sigma_x(Imax + 1 - i);
	alpha_x(i) = alpha_x(Imax + 1 - i);
	%  alpha_x(i) = max_alpha_x;
	RAe_x(i) = (2.0 * EPSILON0 + EM_DELTAT * alpha_x(i)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_x(i) + sigma_x(i))) - 1.0;
	RBe_x(i) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_x(i) + sigma_x(i)));
	RCe_x(i) = EM_DELTAT * (alpha_x(i) - alpha_x(i) * 1.0 - sigma_x(i)) / EPSILON0;
	RDe_x(i) = EM_DELTAT * (alpha_x(i) + sigma_x(i)) / EPSILON0;
	REe_x(i) = 1.0 - RDe_x(i) * RBe_x(i);
	RFe_x(i) = RCe_x(i) - RAe_x(i) * RDe_x(i);
end

% at front boundary (y)
for j = 1:NPML
    sigma_y(j) = max_sigma_y * (((NPML) - j) / (NPML - 1.0))^alpha;
	alpha_y(j) = max_alpha_y * (((NPML) - j) / (NPML - 1.0))^alpha_aa;
	% alpha_y(j) = max_alpha_y;
	RAe_y(j) = (2.0 * EPSILON0 + EM_DELTAT * alpha_y(j)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_y(j) + sigma_y(j))) - 1.0;
	RBe_y(j) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_y(j) + sigma_y(j)));
	RCe_y(j) = EM_DELTAT * (alpha_y(j) - alpha_y(j) * 1.0 - sigma_y(j)) / EPSILON0;
	RDe_y(j) = EM_DELTAT * (alpha_y(j) + sigma_y(j)) / EPSILON0;
	REe_y(j) = 1.0 - RDe_y(j) * RBe_y(j);
	RFe_y(j) = RCe_y(j) - RAe_y(j) * RDe_y(j);
end

% at back boundary (y) 
for j = Jmax - NPML + 1: Jmax
	sigma_y(j) = sigma_y(Jmax + 1 - j);
	alpha_y(j) = alpha_y(Jmax + 1 - j);
	% alpha_y(j) = max_alpha_y;
	RAe_y(j) = (2.0 * EPSILON0 + EM_DELTAT * alpha_y(j)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_y(j) + sigma_y(j))) - 1.0;
	RBe_y(j) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_y(j) + sigma_y(j)));
	RCe_y(j)= EM_DELTAT * (alpha_y(j) - alpha_y(j) * 1.0 - sigma_y(j)) / EPSILON0;
	RDe_y(j) = EM_DELTAT * (alpha_y(j) + sigma_y(j)) / EPSILON0;
	REe_y(j) = 1.0 - RDe_y(j) * RBe_y(j);
	RFe_y(j) = RCe_y(j) - RAe_y(j) * RDe_y(j);
end

% at front boundary (z)
for k = 1:NPML
    sigma_z(k) = max_sigma_z * (((NPML) - k) / (NPML - 1.0))^alpha;
	alpha_z(k) = max_alpha_z * (((NPML) - k) / (NPML - 1.0))^alpha_aa;
	% alpha_z(k) = max_alpha_z;
	RAe_z(k) = (2.0 * EPSILON0 + EM_DELTAT * alpha_z(k)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_z(k) + sigma_z(k))) - 1.0;
	RBe_z(k) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_z(k) + sigma_z(k)));
	RCe_z(k) = EM_DELTAT * (alpha_z(k) - alpha_z(k) * 1.0 - sigma_z(k)) / EPSILON0;
	RDe_z(k) = EM_DELTAT * (alpha_z(k) + sigma_z(k)) / EPSILON0;
	REe_z(k) = 1.0 - RDe_z(k) * RBe_z(k);
	RFe_z(k) = RCe_z(k) - RAe_z(k) * RDe_z(k);
end

% at back boundary (z)
for k = Kmax - NPML + 1: Kmax
	sigma_z(k) = sigma_z(Kmax + 1 - k);
	alpha_z(k) = alpha_z(Kmax + 1 - k);
	% alpha_z(k) = max_alpha_z;
	RAe_z(k) = (2.0 * EPSILON0 + EM_DELTAT * alpha_z(k)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_z(k) + sigma_z(k))) - 1.0;
	RBe_z(k) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_z(k) + sigma_z(k)));
	RCe_z(k) = EM_DELTAT * (alpha_z(k) - alpha_z(k) * 1.0 - sigma_z(k)) / EPSILON0;
	RDe_z(k) = EM_DELTAT * (alpha_z(k) + sigma_z(k)) / EPSILON0;
	REe_z(k) = 1.0 - RDe_z(k) * RBe_z(k);
	RFe_z(k) = RCe_z(k) - RAe_z(k) * RDe_z(k);
end

%%  H  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at front boundary (x)
for i = 1:NPML-1
    sigma_mx(i) = max_sigma_x * (((NPML) - i - 0.5) / (NPML - 1.0))^alpha;
	alpha_mx(i) = max_alpha_x * (((NPML) - i - 0.5) / (NPML - 1.0))^alpha_aa;
	% alpha_x(i) = max_alpha_x;
	RAh_x(i) = (2.0 * EPSILON0 + EM_DELTAT * alpha_mx(i)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_mx(i) + sigma_mx(i))) - 1.0;
	RBh_x(i) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_mx(i) + sigma_mx(i)));
	RCh_x(i) = EM_DELTAT * (alpha_mx(i) - alpha_mx(i) * 1.0 - sigma_mx(i)) / EPSILON0;
	RDh_x(i) = EM_DELTAT * (alpha_mx(i) + sigma_mx(i)) / EPSILON0;
	REh_x(i) = 1.0 - RDh_x(i) * RBh_x(i);
	RFh_x(i) = RCh_x(i) - RAh_x(i) * RDh_x(i);
end

% at back boundary (x)
for i = Imax - NPML + 1 :  Imax - 1
	sigma_mx(i) = sigma_mx(Imax - i);
	alpha_mx(i) = alpha_mx(Imax - i);
	% alpha_x(i) = max_alpha_x;
	RAh_x(i) = (2.0 * EPSILON0 + EM_DELTAT * alpha_mx(i)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_mx(i) + sigma_mx(i))) - 1.0;
	RBh_x(i) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_mx(i) + sigma_mx(i)));
	RCh_x(i) = EM_DELTAT * (alpha_mx(i) - alpha_mx(i) * 1.0 - sigma_mx(i)) / EPSILON0;
	RDh_x(i) = EM_DELTAT * (alpha_mx(i) + sigma_mx(i)) / EPSILON0;
	REh_x(i) = 1.0 - RDh_x(i) * RBh_x(i);
	RFh_x(i) = RCh_x(i) - RAh_x(i) * RDh_x(i);
end

% at front boundary (y)
for j = 1:NPML-1
    sigma_my(j) = max_sigma_y * (((NPML) - j - 0.5) / (NPML - 1.0))^alpha;
	alpha_my(j) = max_alpha_y * (((NPML) - j - 0.5) / (NPML - 1.0))^alpha_aa;
	% alpha_y(j) = max_alpha_y;
	RAh_y(j) = (2.0 * EPSILON0 + EM_DELTAT * alpha_my(j)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_my(j) + sigma_my(j))) - 1.0;
	RBh_y(j) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_my(j) + sigma_my(j)));
	RCh_y(j) = EM_DELTAT * (alpha_my(j) - alpha_my(j) * 1.0 - sigma_my(j)) / EPSILON0;
	RDh_y(j) = EM_DELTAT * (alpha_my(j) + sigma_my(j)) / EPSILON0;
	REh_y(j) = 1.0 - RDh_y(j) * RBh_y(j);
	RFh_y(j) = RCh_y(j) - RAh_y(j) * RDh_y(j);
end

% at back boundary (y)
for j = Jmax - NPML + 1 :  Jmax - 1
	sigma_my(j) = sigma_my(Jmax - j)  ;
	alpha_my(j) = alpha_my(Jmax - j);
	% alpha_y(j) = max_alpha_y;
	RAh_y(j) = (2.0 * EPSILON0 + EM_DELTAT * alpha_my(j)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_my(j) + sigma_my(j))) - 1.0;
	RBh_y(j) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_my(j) + sigma_my(j)));
	RCh_y(j) = EM_DELTAT * (alpha_my(j) - alpha_my(j) * 1.0 - sigma_my(j)) / EPSILON0;
	RDh_y(j) = EM_DELTAT * (alpha_my(j) + sigma_my(j)) / EPSILON0;
	REh_y(j) = 1.0 - RDh_y(j) * RBh_y(j);
	RFh_y(j) = RCh_y(j) - RAh_y(j) * RDh_y(j);
end

% at front boundary (z)
for k = 1:NPML-1
    sigma_mz(k) = max_sigma_z * (((NPML) - k - 0.5) / (NPML - 1.0))^alpha;
	alpha_mz(k) = max_alpha_z * (((NPML) - k - 0.5) / (NPML - 1.0))^alpha_aa;
	% alpha_z(k) = max_alpha_z;
	RAh_z(k) = (2.0 * EPSILON0 + EM_DELTAT * alpha_mz(k)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_mz(k) + sigma_mz(k))) - 1.0;
	RBh_z(k) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_mz(k) + sigma_mz(k)));
	RCh_z(k) = EM_DELTAT * (alpha_mz(k) - alpha_mz(k) * 1.0 - sigma_mz(k)) / EPSILON0;
	RDh_z(k) = EM_DELTAT * (alpha_mz(k) + sigma_mz(k)) / EPSILON0;
	REh_z(k) = 1.0 - RDh_z(k) * RBh_z(k);
	RFh_z(k) = RCh_z(k) - RAh_z(k) * RDh_z(k);
end

% at back boundary (z)
for k = Kmax - NPML + 1 :  Kmax - 1
	sigma_mz(k) = sigma_mz(Kmax - k);
	alpha_mz(k) = alpha_mz(Kmax - k);
% 	alpha_z(k) = max_alpha_z;
	RAh_z(k) = (2.0 * EPSILON0 + EM_DELTAT * alpha_mz(k)) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_mz(k) + sigma_mz(k))) - 1.0;
	RBh_z(k) = (2.0 * EPSILON0) / (2.0 * EPSILON0 + EM_DELTAT * (alpha_mz(k) + sigma_mz(k)));
	RCh_z(k) = EM_DELTAT * (alpha_mz(k) - alpha_mz(k) * 1.0 - sigma_mz(k)) / EPSILON0;
	RDh_z(k) = EM_DELTAT * (alpha_mz(k) + sigma_mz(k)) / EPSILON0;
	REh_z(k) = 1.0 - RDh_z(k) * RBh_z(k);
	RFh_z(k) = RCh_z(k) - RAh_z(k) * RDh_z(k);
end



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%  INITIALIZE VARIABLES 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
%{
sig(:,:,:) = sigM1; 
eps(:,:,:) = epsR*epsO; 

sigk(:,:,:) = sigM1; 
epsk(:,:,:) = epsR*epsO; 
%}

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%  FILL IN UPDATING COEFFICIENTS 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
CA = 1.0; 
CB = (dt/(eps0)); 
 
CP = 1.0; 
CQ = (dt/(mu0)); 
   %{
   for i = 1:Imax 
      for j = 1:Jmax 
         for k = 1:Kmax 
            CA(i,j,k) = (1.0 - sig(i,j,k)*dt / (2.0*eps(i,j,k))) / (1.0 + sig(i,j,k) * dt / (2.0*eps(i,j,k))); 
            CB(i,j,k) = (dt/(eps(i,j,k))) / (1.0 + sig(i,j,k)*dt / (2.0*eps(i,j,k))); 
         end 
      end 
   end 
   
   for i = 1:Imaxk 
      for j = 1:Jmaxk 
         for k = 1:Kmaxk 
            CAk(i,j,k) = (1.0 - sigk(i,j,k)*dt / (2.0*epsk(i,j,k))) / (1.0 + sigk(i,j,k) * dt / (2.0*epsk(i,j,k))); 
            CBk(i,j,k) = (dt/(epsk(i,j,k))) / (1.0 + sigk(i,j,k)*dt / (2.0*epsk(i,j,k))); 
         end 
      end 
   end 
   %}

% Parameters of Debye model for dispaesive media

% Skin
e0_inf_S = 29.9;
es_S = 47.9;
tao_S = 43.6e-12;
deltae_S = es_S - e0_inf_S;
sigma_s_S = 0.540;

CA_S = (eps0 / dt - sigma_s_S / 2.0) / (eps0 / dt + sigma_s_S / 2.0);
CB_S = 1.0 / (eps0 / dt + sigma_s_S / 2.0);

K_S = e0_inf_S;
sigma_S = deltae_S;
Alfa_S = 1.0;

RA_S = ((2.0 * tao_S + dt * Alfa_S) / (2.0 * tao_S * K_S + dt * (Alfa_S * K_S + sigma_S))) - 1.0;
RB_S = (2.0 * tao_S * K_S) / (2.0 * tao_S * K_S + dt * (Alfa_S * K_S + sigma_S));
RC_S = dt * (Alfa_S - Alfa_S * K_S - sigma_S) / (tao_S * K_S);
RD_S = dt * (Alfa_S * K_S + sigma_S) / (tao_S * K_S);
RE_S = 1.0 - RD_S * RB_S;
RF_S = RC_S - RA_S * RD_S;

% Fat
e0_inf_F = 4.0;
es_F = 5.53;
tao_F = 23.6e-12;
deltae_F = es_F - e0_inf_F;
sigma_s_F = 0.037;

CA_F = (eps0 / dt - sigma_s_F / 2.0) / (eps0 / dt + sigma_s_F / 2.0);
CB_F = 1.0 / (eps0 / dt + sigma_s_F / 2.0);

 K_F = e0_inf_F;
sigma_F = deltae_F;
Alfa_F = 1.0;

RA_F = ((2.0 * tao_F + dt * Alfa_F) / (2.0 * tao_F * K_F + dt * (Alfa_F * K_F + sigma_F))) - 1.0;
RB_F = (2.0 * tao_F * K_F) / (2 * tao_F * K_F + dt * (Alfa_F * K_F + sigma_F));
RC_F = dt * (Alfa_F - Alfa_F * K_F - sigma_F) / (tao_F * K_F);
RD_F = dt * (Alfa_F * K_F + sigma_F) / (tao_F * K_F);
RE_F = 1.0 - RD_F * RB_F;
RF_F = RC_F - RA_F * RD_F;

% Bone
e0_inf_M = 7.36;
es_M = 14.2;
tao_M = 34.1e-12;
deltae_M = es_M - e0_inf_M;
sigma_s_M = 0.104;

CA_M = (eps0 / dt - sigma_s_M / 2.0) / (eps0 / dt + sigma_s_M / 2.0);
CB_M = 1.0 / (eps0 / dt + sigma_s_M / 2.0);

K_M = e0_inf_M;
sigma_M = deltae_M;
Alfa_M = 1.0;

RA_M = ((2.0 * tao_M + dt * Alfa_M) / (2.0 * tao_M * K_M + dt * (Alfa_M * K_M + sigma_M))) - 1.0;
RB_M = (2.0 * tao_M * K_M) / (2.0 * tao_M * K_M + dt * (Alfa_M * K_M + sigma_M));
RC_M = dt * (Alfa_M - Alfa_M * K_M - sigma_M) / (tao_M * K_M);
RD_M = dt * (Alfa_M * K_M + sigma_M) / (tao_M * K_M);
RE_M = 1.0 - RD_M * RB_M;
RF_M = RC_M - RA_M * RD_M;
                                                                      
IC = round((Imax + 1) / 2.0) + 1;
JC = round(1.2*(Jmax + 1) / 2.0) + 0;
KC = round((Kmax + 1) / 2.0) + 1;

R = 96;
R1 = 35;
HL = KC - 25;
HP = KC + 25;
zaj = JC;
zbj = zaj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
		    CA_Tx(i,j,k) = CA;
		  	CB_Tx(i,j,k) = CB;
            
            if (i >= IC - R1 && i <= IC + R1 && j >= zbj && j <= zbj && k >= HL && k <= HP)
				CA_Tx(i,j,k) = CA_S;
				CB_Tx(i,j,k) = CB_S;
            elseif (i >= IC - R1 && i <= IC + R1 && j >= zbj + R && j <= zbj + R && k >= HL && k <= HP)
				CA_Tx(i,j,k) = CA_S;
				CB_Tx(i,j,k) = CB_S;
            elseif (i >= IC - R1 && i <= IC + R1 && j >= zbj + 1 && j <= zbj + 2 && k >= HL && k <= HP)   % 左
				CA_Tx(i,j,k) = CA_F;
				CB_Tx(i,j,k) = CB_F;
             elseif (i >= IC - R1 && i <= IC + R1 && j >= (zbj + R - 1) - 1 && j <= (zbj + R - 1) && k >= HL && k <= HP)  % 右		
                CA_Tx(i,j,k) = CA_F;
				CB_Tx(i,j,k) = CB_F;
            elseif (i >= IC - R1 && i <= IC - R1 + 2 && j >= zbj + 3 && j <= (zbj + R - 1) - 2 && k >= HL && k <= HP)  % 上
				CA_Tx(i,j,k) = CA_F;
				CB_Tx(i,j,k) = CB_F;
            elseif (i >= IC + R1 - 2 && i <= IC + R1 && j >= zbj + 3 && j <= (zbj + R - 1) - 2 && k >= HL && k <= HP)  % 下
				CA_Tx(i,j,k) = CA_F;
				CB_Tx(i,j,k) = CB_F;
            elseif (i >= IC - R1 +1 + 2 && i <= IC + R1 - 1 -2 && j >= zbj + 3 && j <= (zbj + R) - 3 && k >= HL && k <= HP)  % 内部
				CA_Tx(i,j,k) = CA_M;
				CB_Tx(i,j,k) = CB_M;
            end
        end
    end
end

for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
			 CA_Ty(i,j,k) = CA;
			 CB_Ty(i,j,k) = CB;
             
             if (i >= IC - R1 && i <= IC - R1 && j >= JC && j <= JC + R -1 && k >= HL && k <= HP)
				CA_Ty(i,j,k) = CA_S;
				CB_Ty(i,j,k) = CB_S;
             elseif (i >= IC + R1 + 1 && i <= IC + R1 + 1 && j >= JC && j <= JC + R - 1 && k >= HL && k <= HP)
				CA_Ty(i,j,k) = CA_S;
				CB_Ty(i,j,k) = CB_S;
             elseif (i >= IC - R1 + 1 && i <= IC - R1 + 2 && j >= JC && j <= JC + R - 1 && k >= HL && k <= HP)   % 上
				CA_Ty(i,j,k) = CA_F;
				CB_Ty(i,j,k) = CB_F;
             elseif (i >= (IC + R1 + 1) - 2 && i <= (IC + R1 + 1) - 1 && j >= JC && j <= JC + R - 1 && k >= HL && k <= HP)  % 下
				CA_Ty(i,j,k) = CA_F;
				CB_Ty(i,j,k) = CB_F;
             elseif (i >= IC - R1 + 3 && i <= (IC + R1 + 1) - 3 && j >= JC && j <= JC + 2 && k >= HL && k <= HP)  % 左
				CA_Ty(i,j,k) = CA_F;
				CB_Ty(i,j,k) = CB_F;
             elseif (i >= IC - R1 + 3 && i <= (IC + R1 + 1) - 3 && j >= JC + R - 3 && j <= JC + R - 1 && k >= HL && k <= HP)  % 右
				CA_Ty(i,j,k) = CA_F;
				CB_Ty(i,j,k) = CB_F;
             elseif (i >= IC - R1 + 3 && i <= (IC + R1 + 1) - 3 && j >= JC + 3 && j <= JC + R - 4 && k >= HL && k <= HP)  % 内部
				CA_Ty(i,j,k) = CA_M;
				CB_Ty(i,j,k) = CB_M;
             end
        end
    end
end

for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
			 CA_Tz(i,j,k) = CA;
			 CB_Tz(i,j,k) = CB;
             
             if (i >= IC - R1 && i <= IC + R1 + 1 && j >= JC && j <= JC && k >= HL && k <= HP - 1)  % 左边
				CA_Tz(i,j,k) = CA_S;
				CB_Tz(i,j,k) = CB_S;
             elseif (i >= IC - R1 && i <= IC + R1 + 1 && j >= zbj + R && j <= zbj + R  && k >= HL && k <= HP - 1) % 右边
				CA_Tz(i,j,k) = CA_S;
				CB_Tz(i,j,k) = CB_S;
             elseif (i >= IC - R1 && i <= IC - R1 && j >= JC + 1 && j <= JC + R - 1 && k >= HL && k <= HP - 1)   % 上边
				CA_Tz(i,j,k) = CA_S;
				CB_Tz(i,j,k) = CB_S;
             elseif (i >= (IC + R1 + 1) && i <= (IC + R1 + 1) && j >= JC + 1 && j <= JC + R - 1 && k >= HL && k <= HP - 1)  % 下边
				CA_Tz(i,j,k) = CA_S;
				CB_Tz(i,j,k) = CB_S;
             elseif (i >= IC - R1 + 1 && i <= (IC + R1 + 1) - 1 && j >= JC + 1 && j <= JC + 2 && k >= HL && k <= HP - 1)  % 左边
				CA_Tz(i,j,k) = CA_F;
				CB_Tz(i,j,k) = CB_F;
             elseif (i >= IC - R1 + 1 && i <= (IC + R1 + 1) - 1 && j >= JC + R - 2 && j <= JC + R - 1 && k >= HL && k <= HP - 1)  % 右边
				CA_Tz(i,j,k) = CA_F;
				CB_Tz(i,j,k) = CB_F;
             elseif (i >= IC - R1 + 1 && i <= (IC - R1 + 1) + 1 && j >= JC + 3 && j <= JC + R - 3 && k >= HL && k <= HP - 1)  % 上边
				CA_Tz(i,j,k) = CA_F;
				CB_Tz(i,j,k)= CB_F;
             elseif (i >= IC + R1 - 1 && i <= (IC + R1 - 1) + 1 && j >= JC + 3 && j <= JC + R - 3 && k >= HL && k <= HP - 1)  % 下边
				CA_Tz(i,j,k) = CA_F;
				CB_Tz(i,j,k) = CB_F;
             elseif (i >= IC - R1 + 3 && i <= (IC + R1 + 1) - 3 && j >= JC + 3 && j <= JC + R - 3 && k >= HL && k <= HP - 1)  % 内部
				CA_Tz(i,j,k) = CA_M;
				CB_Tz(i,j,k) = CB_M;
             end
        end
    end
end
	
for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
            
            CB_index_X(i,j,k) = 0.0;
            
			RA_index_X_y(i,j,k) = RAe_y(j);
			RA_index_X_z(i,j,k) = RAe_z(k);

			RB_index_X_y(i,j,k) = RBe_y(j);
			RB_index_X_z(i,j,k) = RBe_z(k);

			RE_index_X_y(i,j,k) = REe_y(j);
			RE_index_X_z(i,j,k) = REe_z(k);

			RF_index_X_y(i,j,k) = RFe_y(j);
			RF_index_X_z(i,j,k) = RFe_z(k);

			RA_index_X(i,j,k) = 0.0;
			RB_index_X(i,j,k) = 0.0;
			RE_index_X(i,j,k) = 0.0;
			RF_index_X(i,j,k) = 0.0;
			sigma_s_index_X(i,j,k) = 0.0;
            
            if (j >= 1 && j < NPML + 1)
				CB_index_X(i,j,k) = CB;
                elseif (j >= Jmax - NPML + 1 && j < Jmax + 1)
				CB_index_X(i,j,k) = CB;
                elseif (k >= 1 && k < NPML + 1)
				CB_index_X(i,j,k) = CB;
                elseif (k >= Kmax - NPML + 1 && k < Kmax + 1)
				CB_index_X(i,j,k) = CB;
            end
            
            if (i >= IC - R1 && i <= IC + R1 && j >= zbj && j <= zbj && k >= HL && k <= HP)
                
				CB_index_X(i,j,k) = CB_S;

				RA_index_X_y(i,j,k) = RA_S;
				RA_index_X_z(i,j,k) = RA_S;

				RB_index_X_y(i,j,k) = RB_S;
				RB_index_X_z(i,j,k) = RB_S;

				RE_index_X_y(i,j,k) = RE_S;
				RE_index_X_z(i,j,k) = RE_S;

				RF_index_X_y(i,j,k) = RF_S;
				RF_index_X_z(i,j,k) = RF_S;

				RA_index_X(i,j,k) = RA_S;
				RB_index_X(i,j,k) = RB_S;
				RE_index_X(i,j,k) = RE_S;
				RF_index_X(i,j,k) = RF_S;
				sigma_s_index_X(i,j,k) = sigma_s_S;
            elseif (i >= IC - R1 && i <= IC + R1 && j >= zbj + R && j <= zbj + R && k >= HL && k <= HP)
				CB_index_X(i,j,k) = CB_S;

				RA_index_X_y(i,j,k) = RA_S;
				RA_index_X_z(i,j,k) = RA_S;

				RB_index_X_y(i,j,k) = RB_S;
				RB_index_X_z(i,j,k) = RB_S;

				RE_index_X_y(i,j,k) = RE_S;
				RE_index_X_z(i,j,k) = RE_S;

				RF_index_X_y(i,j,k) = RF_S;
				RF_index_X_z(i,j,k) = RF_S;

				RA_index_X(i,j,k) = RA_S;
				RB_index_X(i,j,k) = RB_S;
				RE_index_X(i,j,k) = RE_S;
				RF_index_X(i,j,k) = RF_S;
				sigma_s_index_X(i,j,k) = sigma_s_S;
            elseif (i >= IC - R1 && i <= IC + R1 && j >= zbj + 1 && j <= zbj + 2 && k >= HL && k <= HP)   % 左
				CB_index_X(i,j,k) = CB_F;

				RA_index_X_y(i,j,k) = RA_F;
				RA_index_X_z(i,j,k) = RA_F;

				RB_index_X_y(i,j,k) = RB_F;
				RB_index_X_z(i,j,k) = RB_F;

				RE_index_X_y(i,j,k) = RE_F;
				RE_index_X_z(i,j,k) = RE_F;

				RF_index_X_y(i,j,k) = RF_F;
				RF_index_X_z(i,j,k) = RF_F;

				RA_index_X(i,j,k) = RA_F;
				RB_index_X(i,j,k) = RB_F;
				RE_index_X(i,j,k) = RE_F;
				RF_index_X(i,j,k) = RF_F;
				sigma_s_index_X(i,j,k) = sigma_s_F;
            elseif (i >= IC - R1 && i <= IC + R1 && j >= (zbj + R - 1) - 1 && j <= (zbj + R - 1) && k >= HL && k <= HP)  % 右
					
                CB_index_X(i,j,k) = CB_F;

				RA_index_X_y(i,j,k) = RA_F;
				RA_index_X_z(i,j,k) = RA_F;

				RB_index_X_y(i,j,k) = RB_F;
				RB_index_X_z(i,j,k) = RB_F;

				RE_index_X_y(i,j,k) = RE_F;
				RE_index_X_z(i,j,k) = RE_F;

				RF_index_X_y(i,j,k) = RF_F;
				RF_index_X_z(i,j,k) = RF_F;

				RA_index_X(i,j,k) = RA_F;
				RB_index_X(i,j,k) = RB_F;
				RE_index_X(i,j,k) = RE_F;
				RF_index_X(i,j,k) = RF_F;
				sigma_s_index_X(i,j,k) = sigma_s_F;
            elseif (i >= IC - R1 && i <= IC - R1 + 2 && j >= zbj + 3 && j <= (zbj + R - 1) - 2 && k >= HL && k <= HP)  % 上
                
					CB_index_X(i,j,k) = CB_F;

					RA_index_X_y(i,j,k) = RA_F;
					RA_index_X_z(i,j,k) = RA_F;

					RB_index_X_y(i,j,k) = RB_F;
					RB_index_X_z(i,j,k) = RB_F;

					RE_index_X_y(i,j,k) = RE_F;
					RE_index_X_z(i,j,k) = RE_F;

					RF_index_X_y(i,j,k) = RF_F;
					RF_index_X_z(i,j,k) = RF_F;

					RA_index_X(i,j,k) = RA_F;
					RB_index_X(i,j,k) = RB_F;
					RE_index_X(i,j,k) = RE_F;
					RF_index_X(i,j,k) = RF_F;
					sigma_s_index_X(i,j,k) = sigma_s_F;
            elseif (i >= IC + R1 - 2 && i <= IC + R1 && j >= zbj + 3 && j <= (zbj + R - 1) - 2 && k >= HL && k <= HP)  % 下
					CB_index_X(i,j,k) = CB_F;

					RA_index_X_y(i,j,k) = RA_F;
					RA_index_X_z(i,j,k) = RA_F;

					RB_index_X_y(i,j,k) = RB_F;
					RB_index_X_z(i,j,k) = RB_F;

					RE_index_X_y(i,j,k) = RE_F;
					RE_index_X_z(i,j,k) = RE_F;

					RF_index_X_y(i,j,k) = RF_F;
					RF_index_X_z(i,j,k) = RF_F;

					RA_index_X(i,j,k) = RA_F;
					RB_index_X(i,j,k) = RB_F;
					RE_index_X(i,j,k) = RE_F;
					RF_index_X(i,j,k) = RF_F;
					sigma_s_index_X(i,j,k) = sigma_s_F;
            elseif (i >= IC - R1 +1 + 2 && i <= IC + R1 - 1 -2 && j >= zbj + 3 && j <= (zbj + R) - 3 && k >= HL && k <= HP)  % 内部
					CB_index_X(i,j,k) = CB_M;

					RA_index_X_y(i,j,k) = RA_M;
					RA_index_X_z(i,j,k) = RA_M;

					RB_index_X_y(i,j,k) = RB_M;
					RB_index_X_z(i,j,k) = RB_M;

					RE_index_X_y(i,j,k) = RE_M;
					RE_index_X_z(i,j,k) = RE_M;

					RF_index_X_y(i,j,k) = RF_M;
					RF_index_X_z(i,j,k) = RF_M;

					RA_index_X(i,j,k) = RA_M;
					RB_index_X(i,j,k) = RB_M;
					RE_index_X(i,j,k) = RE_M;
					RF_index_X(i,j,k) = RF_M;
					sigma_s_index_X(i,j,k) = sigma_s_M;
            end
            
        end
    end
end

for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
            
            CB_index_Y(i,j,k) = 0.0;
            
			RA_index_Y_x(i,j,k) = RAe_x(i);
			RA_index_Y_z(i,j,k) = RAe_z(k);

			RB_index_Y_x(i,j,k) = RBe_x(i);
			RB_index_Y_z(i,j,k) = RBe_z(k);

			RE_index_Y_x(i,j,k) = REe_x(i);
			RE_index_Y_z(i,j,k) = REe_z(k);

			RF_index_Y_x(i,j,k) = RFe_x(i);
			RF_index_Y_z(i,j,k) = RFe_z(k);

			RA_index_Y(i,j,k) = 0.0;
			RB_index_Y(i,j,k) = 0.0;
			RE_index_Y(i,j,k) = 0.0;
			RF_index_Y(i,j,k) = 0.0;
			sigma_s_index_Y(i,j,k) = 0.0;
            
             if (i >= 1 && i < NPML + 1)
				CB_index_Y(i,j,k) = CB;
            elseif (i >= Imax - NPML + 1 && i < Imax + 1)
				CB_index_Y(i,j,k) = CB;
            elseif (k >= 1 && k < NPML + 1)
				CB_index_Y(i,j,k) = CB;
            elseif (k >= Kmax - NPML + 1 && k < Kmax + 1)
				CB_index_Y(i,j,k) = CB;
             end
            
            if (i >= IC - R1 && i <= IC - R1 && j >= JC && j <= JC + R -1 && k >= HL && k <= HP)
				CB_index_Y(i,j,k) = CB_S;

				RA_index_Y_x(i,j,k) = RA_S;
				RA_index_Y_z(i,j,k) = RA_S;

				RB_index_Y_x(i,j,k) = RB_S;
				RB_index_Y_z(i,j,k) = RB_S;

				RE_index_Y_x(i,j,k) = RE_S;
				RE_index_Y_z(i,j,k) = RE_S;

				RF_index_Y_x(i,j,k) = RF_S;
				RF_index_Y_z(i,j,k) = RF_S;

				RA_index_Y(i,j,k) = RA_S;
				RB_index_Y(i,j,k) = RB_S;
				RE_index_Y(i,j,k) = RE_S;
				RF_index_Y(i,j,k) = RF_S;
				sigma_s_index_Y(i,j,k) = sigma_s_S;
            elseif (i >= IC + R1 + 1 && i <= IC + R1 + 1 && j >= JC && j <= JC + R - 1 && k >= HL && k <= HP)
					CB_index_Y(i,j,k) = CB_S;

					RA_index_Y_x(i,j,k) = RA_S;
					RA_index_Y_z(i,j,k) = RA_S;

					RB_index_Y_x(i,j,k) = RB_S;
					RB_index_Y_z(i,j,k) = RB_S;

					RE_index_Y_x(i,j,k) = RE_S;
					RE_index_Y_z(i,j,k) = RE_S;

					RF_index_Y_x(i,j,k) = RF_S;
					RF_index_Y_z(i,j,k) = RF_S;

					RA_index_Y(i,j,k) = RA_S;
					RB_index_Y(i,j,k) = RB_S;
					RE_index_Y(i,j,k) = RE_S;
					RF_index_Y(i,j,k) = RF_S;
					sigma_s_index_Y(i,j,k) = sigma_s_S;
            elseif (i >= IC - R1 + 1 && i <= IC - R1 + 2 && j >= JC && j <= JC + R - 1 && k >= HL && k <= HP)   % 上
					CB_index_Y(i,j,k) = CB_F;

					RA_index_Y_x(i,j,k) = RA_F;
					RA_index_Y_z(i,j,k) = RA_F;

					RB_index_Y_x(i,j,k) = RB_F;
					RB_index_Y_z(i,j,k) = RB_F;

					RE_index_Y_x(i,j,k) = RE_F;
					RE_index_Y_z(i,j,k) = RE_F;

					RF_index_Y_x(i,j,k) = RF_F;
					RF_index_Y_z(i,j,k) = RF_F;

					RA_index_Y(i,j,k) = RA_F;
					RB_index_Y(i,j,k) = RB_F;
					RE_index_Y(i,j,k) = RE_F;
					RF_index_Y(i,j,k) = RF_F;
					sigma_s_index_Y(i,j,k) = sigma_s_F;
            elseif (i >= (IC + R1 + 1) - 2 && i <= (IC + R1 + 1) - 1 && j >= JC && j <= JC + R - 1 && k >= HL && k <= HP)  % 下
					CB_index_Y(i,j,k) = CB_F;

					RA_index_Y_x(i,j,k) = RA_F;
					RA_index_Y_z(i,j,k) = RA_F;

					RB_index_Y_x(i,j,k) = RB_F;
					RB_index_Y_z(i,j,k) = RB_F;

					RE_index_Y_x(i,j,k) = RE_F;
					RE_index_Y_z(i,j,k) = RE_F;

					RF_index_Y_x(i,j,k) = RF_F;
					RF_index_Y_z(i,j,k) = RF_F;

					RA_index_Y(i,j,k) = RA_F;
					RB_index_Y(i,j,k) = RB_F;
					RE_index_Y(i,j,k) = RE_F;
					RF_index_Y(i,j,k) = RF_F;
					sigma_s_index_Y(i,j,k) = sigma_s_F;
            elseif (i >= IC - R1 + 3 && i <= (IC + R1 + 1) - 3 && j >= JC && j <= JC + 2 && k >= HL && k <= HP)  % 左
					CB_index_Y(i,j,k) = CB_F;

					RA_index_Y_x(i,j,k) = RA_F;
					RA_index_Y_z(i,j,k) = RA_F;

					RB_index_Y_x(i,j,k) = RB_F;
					RB_index_Y_z(i,j,k) = RB_F;

					RE_index_Y_x(i,j,k) = RE_F;
					RE_index_Y_z(i,j,k) = RE_F;

					RF_index_Y_x(i,j,k) = RF_F;
					RF_index_Y_z(i,j,k) = RF_F;

					RA_index_Y(i,j,k) = RA_F;
					RB_index_Y(i,j,k) = RB_F;
					RE_index_Y(i,j,k) = RE_F;
					RF_index_Y(i,j,k) = RF_F;
					sigma_s_index_Y(i,j,k) = sigma_s_F;
            elseif (i >= IC - R1 + 3 && i <= (IC + R1 + 1) - 3 && j >= JC + R - 3 && j <= JC + R - 1 && k >= HL && k <= HP)  % 右
					CB_index_Y(i,j,k) = CB_F;

					RA_index_Y_x(i,j,k) = RA_F;
					RA_index_Y_z(i,j,k) = RA_F;

					RB_index_Y_x(i,j,k) = RB_F;
					RB_index_Y_z(i,j,k) = RB_F;

					RE_index_Y_x(i,j,k) = RE_F;
					RE_index_Y_z(i,j,k) = RE_F;

					RF_index_Y_x(i,j,k) = RF_F;
					RF_index_Y_z(i,j,k) = RF_F;

					RA_index_Y(i,j,k) = RA_F;
					RB_index_Y(i,j,k) = RB_F;
					RE_index_Y(i,j,k) = RE_F;
					RF_index_Y(i,j,k) = RF_F;
					sigma_s_index_Y(i,j,k) = sigma_s_F;
            elseif (i >= IC - R1 + 3 && i <= (IC + R1 + 1) - 3 && j >= JC + 3 && j <= JC + R - 4 && k >= HL && k <= HP)  % 内部
					CB_index_Y(i,j,k) = CB_M;

					RA_index_Y_x(i,j,k) = RA_M;
					RA_index_Y_z(i,j,k) = RA_M;

					RB_index_Y_x(i,j,k) = RB_M;
					RB_index_Y_z(i,j,k) = RB_M;

					RE_index_Y_x(i,j,k) = RE_M;
					RE_index_Y_z(i,j,k) = RE_M;

					RF_index_Y_x(i,j,k) = RF_M;
					RF_index_Y_z(i,j,k) = RF_M;

					RA_index_Y(i,j,k) = RA_M;
					RB_index_Y(i,j,k) = RB_M;
					RE_index_Y(i,j,k) = RE_M;
					RF_index_Y(i,j,k) = RF_M;
					sigma_s_index_Y(i,j,k) = sigma_s_M;
            end
            
        end
    end
end

for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1

            CB_index_Z(i,j,k) = 0.0;
            
			RA_index_Z_x(i,j,k) = RAe_x(i);
			RA_index_Z_y(i,j,k) = RAe_y(j);

			RB_index_Z_x(i,j,k) = RBe_x(i);
			RB_index_Z_y(i,j,k) = RBe_y(j);

			RE_index_Z_x(i,j,k) = REe_x(i);
			RE_index_Z_y(i,j,k) = REe_y(j);

			RF_index_Z_x(i,j,k) = RFe_x(i);
			RF_index_Z_y(i,j,k) = RFe_y(j);

			RA_index_Z(i,j,k) = 0.0;
			RB_index_Z(i,j,k) = 0.0;
			RE_index_Z(i,j,k) = 0.0;
			RF_index_Z(i,j,k) = 0.0;
			sigma_s_index_Z(i,j,k) = 0.0;
            
              if (i >= 1 && i < NPML+1)
			     CB_index_Z(i,j,k) = CB;
            elseif (i >= Imax - NPML + 1 && i < Imax + 1)
		         CB_index_Z(i,j,k) = CB;
            elseif (j >= 1 && j < NPML + 1)
			     CB_index_Z(i,j,k) = CB;
            elseif (j >= Jmax - NPML + 1 && j < Jmax + 1)
			     CB_index_Z(i,j,k) = CB;
             end 
            
            if (i >= IC - R1 && i <= IC + R1 + 1 && j >= JC && j <= JC && k >= HL && k <= HP - 1)  % 左边
			        CB_index_Z(i,j,k) = CB_S;

			        RA_index_Z_x(i,j,k) = RA_S;
			        RA_index_Z_y(i,j,k) = RA_S;
 
			        RB_index_Z_x(i,j,k) = RB_S;
			        RB_index_Z_y(i,j,k) = RB_S;

			        RE_index_Z_x(i,j,k) = RE_S;
			        RE_index_Z_y(i,j,k) = RE_S;

			        RF_index_Z_x(i,j,k) = RF_S;
			        RF_index_Z_y(i,j,k) = RF_S;
                    
                    RA_index_Z(i,j,k) = RA_S;
                    RB_index_Z(i,j,k) = RB_S;
                    RE_index_Z(i,j,k) = RE_S;
                    RF_index_Z(i,j,k) = RF_S;
		            sigma_s_index_Z(i,j,k) = sigma_s_S;
            elseif (i >= IC - R1 && i <= IC + R1 + 1 && j >= zbj + R && j <= zbj + R  && k >= HL && k <= HP - 1) % 右边
			        CB_index_Z(i,j,k) = CB_S;

		            RA_index_Z_x(i,j,k) = RA_S;
			        RA_index_Z_y(i,j,k) = RA_S;

			        RB_index_Z_x(i,j,k) = RB_S;
			        RB_index_Z_y(i,j,k) = RB_S;
 
			        RE_index_Z_x(i,j,k) = RE_S;
			        RE_index_Z_y(i,j,k) = RE_S;

			        RF_index_Z_x(i,j,k) = RF_S;
			        RF_index_Z_y(i,j,k) = RF_S;

			        RA_index_Z(i,j,k) = RA_S;
			        RB_index_Z(i,j,k) = RB_S;
			        RE_index_Z(i,j,k) = RE_S;
			        RF_index_Z(i,j,k) = RF_S;
			        sigma_s_index_Z(i,j,k) = sigma_s_S;
            elseif (i >= IC - R1 && i <= IC - R1 && j >= JC + 1 && j <= JC + R - 1 && k >= HL && k <= HP - 1)   % 上边
			        CB_index_Z(i,j,k) = CB_S;

			        RA_index_Z_x(i,j,k) = RA_S;
			        RA_index_Z_y(i,j,k) = RA_S;

			        RB_index_Z_x(i,j,k) = RB_S;
			        RB_index_Z_y(i,j,k) = RB_S;

			        RE_index_Z_x(i,j,k) = RE_S;
			        RE_index_Z_y(i,j,k) = RE_S;

			        RF_index_Z_x(i,j,k) = RF_S;
			        RF_index_Z_y(i,j,k) = RF_S;

			        RA_index_Z(i,j,k) = RA_S;
			        RB_index_Z(i,j,k) = RB_S;
			        RE_index_Z(i,j,k) = RE_S;
			        RF_index_Z(i,j,k) = RF_S;
			        sigma_s_index_Z(i,j,k) = sigma_s_S;
            elseif (i >= (IC + R1 + 1) && i <= (IC + R1 + 1) && j >= JC + 1 && j <= JC + R - 1 && k >= HL && k <= HP - 1)  % 下边
			        CB_index_Z(i,j,k) = CB_S;

			        RA_index_Z_x(i,j,k) = RA_S;
			        RA_index_Z_y(i,j,k) = RA_S;

			        RB_index_Z_x(i,j,k) = RB_S;
			        RB_index_Z_y(i,j,k) = RB_S;

			        RE_index_Z_x(i,j,k) = RE_S;
			        RE_index_Z_y(i,j,k) = RE_S;

			        RF_index_Z_x(i,j,k) = RF_S;
		      	    RF_index_Z_y(i,j,k) = RF_S;

			        RA_index_Z(i,j,k) = RA_S;
			        RB_index_Z(i,j,k) = RB_S;
			        RE_index_Z(i,j,k) = RE_S;
			        RF_index_Z(i,j,k) = RF_S;
			        sigma_s_index_Z(i,j,k) = sigma_s_S;
            elseif (i >= IC - R1 + 1 && i <= (IC + R1 + 1) - 1 && j >= JC + 1 && j <= JC + 2 && k >= HL && k <= HP - 1)  % 左边
			        CB_index_Z(i,j,k) = CB_F;

			        RA_index_Z_x(i,j,k) = RA_F;
			        RA_index_Z_y(i,j,k) = RA_F;

			        RB_index_Z_x(i,j,k) = RB_F;
			        RB_index_Z_y(i,j,k) = RB_F;

		  	        RE_index_Z_x(i,j,k) = RE_F;
			        RE_index_Z_y(i,j,k) = RE_F;

			        RF_index_Z_x(i,j,k) = RF_F;
			        RF_index_Z_y(i,j,k) = RF_F;

			        RA_index_Z(i,j,k) = RA_F;
			        RB_index_Z(i,j,k) = RB_F;
			        RE_index_Z(i,j,k) = RE_F;
			        RF_index_Z(i,j,k) = RF_F;
			        sigma_s_index_Z(i,j,k) = sigma_s_F;
            elseif (i >= IC - R1 + 1 && i <= (IC + R1 + 1) - 1 && j >= JC + R - 2 && j <= JC + R - 1 && k >= HL && k <= HP - 1)  % 右边
			        CB_index_Z(i,j,k) = CB_F;

			        RA_index_Z_x(i,j,k) = RA_F;
			        RA_index_Z_y(i,j,k) = RA_F;

			        RB_index_Z_x(i,j,k) = RB_F;
			        RB_index_Z_y(i,j,k) = RB_F;

			        RE_index_Z_x(i,j,k) = RE_F;
			        RE_index_Z_y(i,j,k) = RE_F;

			        RF_index_Z_x(i,j,k) = RF_F;
			        RF_index_Z_y(i,j,k) = RF_F;

			        RA_index_Z(i,j,k) = RA_F;
			        RB_index_Z(i,j,k) = RB_F;
			        RE_index_Z(i,j,k) = RE_F;
			        RF_index_Z(i,j,k) = RF_F;
			        sigma_s_index_Z(i,j,k) = sigma_s_F;
            elseif (i >= IC - R1 + 1 && i <= (IC - R1 + 1) + 1 && j >= JC + 3 && j <= JC + R - 3 && k >= HL && k <= HP - 1)  % 上边
			        CB_index_Z(i,j,k) = CB_F;

			        RA_index_Z_x(i,j,k) = RA_F;
			        RA_index_Z_y(i,j,k) = RA_F;

			        RB_index_Z_x(i,j,k) = RB_F;
			        RB_index_Z_y(i,j,k) = RB_F;

			        RE_index_Z_x(i,j,k) = RE_F;
			        RE_index_Z_y(i,j,k) = RE_F;

			        RF_index_Z_x(i,j,k) = RF_F;
			        RF_index_Z_y(i,j,k) = RF_F;

			        RA_index_Z(i,j,k) = RA_F;
			        RB_index_Z(i,j,k) = RB_F;
			        RE_index_Z(i,j,k) = RE_F;
			        RF_index_Z(i,j,k) = RF_F;
			        sigma_s_index_Z(i,j,k) = sigma_s_F;
            elseif (i >= IC + R1 - 1 && i <= (IC + R1 - 1) + 1 && j >= JC + 3 && j <= JC + R - 3 && k >= HL && k <= HP - 1)  % 下边
			        CB_index_Z(i,j,k) = CB_F;

			        RA_index_Z_x(i,j,k) = RA_F;
			        RA_index_Z_y(i,j,k) = RA_F;

			        RB_index_Z_x(i,j,k) = RB_F;
			        RB_index_Z_y(i,j,k) = RB_F;

			        RE_index_Z_x(i,j,k) = RE_F;
			        RE_index_Z_y(i,j,k) = RE_F;

			        RF_index_Z_x(i,j,k) = RF_F;
		            RF_index_Z_y(i,j,k) = RF_F;

			        RA_index_Z(i,j,k) = RA_F;
			        RB_index_Z(i,j,k) = RB_F;
			        RE_index_Z(i,j,k) = RE_F;
			        RF_index_Z(i,j,k) = RF_F;
			        sigma_s_index_Z(i,j,k) = sigma_s_F;
            elseif (i >= IC - R1 + 3 && i <= (IC + R1 + 1) - 3 && j >= JC + 3 && j <= JC + R - 3 && k >= HL && k <= HP - 1)  % 内部
			        CB_index_Z(i,j,k) = CB_M;

			        RA_index_Z_x(i,j,k) = RA_M;
			        RA_index_Z_y(i,j,k) = RA_M;

			        RB_index_Z_x(i,j,k) = RB_M;
			        RB_index_Z_y(i,j,k) = RB_M;

			        RE_index_Z_x(i,j,k) = RE_M;
			        RE_index_Z_y(i,j,k) = RE_M;

			        RF_index_Z_x(i,j,k) = RF_M;
			        RF_index_Z_y(i,j,k) = RF_M;

			        RA_index_Z(i,j,k) = RA_M;
			        RB_index_Z(i,j,k) = RB_M;
			        RE_index_Z(i,j,k) = RE_M;
			        RF_index_Z(i,j,k) = RF_M;
			       sigma_s_index_Z(i,j,k) = sigma_s_M;
            end
     
        end
    end
end  

CB_RA_X_y = zeros(Imax+1, Jmax+1, Kmax+1);
CB_RA_X_z = zeros(Imax+1, Jmax+1, Kmax+1);
                
CB_RB_X_y = zeros(Imax+1, Jmax+1, Kmax+1);
CB_RB_X_z = zeros(Imax+1, Jmax+1, Kmax+1);
             
CB_RA_sigma_s_X_X = zeros(Imax+1, Jmax+1, Kmax+1);
CB_RB_sigma_s_X_X = zeros(Imax+1, Jmax+1, Kmax+1);

 for i=1:Imax+1
	 for j=1:Jmax+1
	     for k=1:Kmax+1
                
             CB_RA_X_y(i,j,k)= CB_index_X(i,j,k) * RA_index_X_y(i,j,k);  
             CB_RA_X_z(i,j,k)= CB_index_X(i,j,k) * RA_index_X_z(i,j,k);  
                
             CB_RB_X_y(i,j,k)= CB_index_X(i,j,k) *  RB_index_X_y(i,j,k);  
             CB_RB_X_z(i,j,k)= CB_index_X(i,j,k) *  RB_index_X_z(i,j,k);  
             
             CB_RA_sigma_s_X_X(i,j,k)= sigma_s_index_X(i,j,k) * CB_index_X(i,j,k) * RA_index_X(i,j,k);  
             CB_RB_sigma_s_X_X(i,j,k)= sigma_s_index_X(i,j,k) * CB_index_X(i,j,k) * RB_index_X(i,j,k);
             
         end
     end
 end

CB_RA_Y_x = zeros(Imax+1, Jmax+1, Kmax+1);
CB_RA_Y_z = zeros(Imax+1, Jmax+1, Kmax+1);
                
CB_RB_Y_x = zeros(Imax+1, Jmax+1, Kmax+1);
CB_RB_Y_z = zeros(Imax+1, Jmax+1, Kmax+1);
             
CB_RA_sigma_s_Y_Y = zeros(Imax+1, Jmax+1, Kmax+1);
CB_RB_sigma_s_Y_Y = zeros(Imax+1, Jmax+1, Kmax+1);

 for i=1:Imax+1
	 for j=1:Jmax+1
	     for k=1:Kmax+1
                
             CB_RA_Y_x(i,j,k)= CB_index_Y(i,j,k) * RA_index_Y_x(i,j,k);  
             CB_RA_Y_z(i,j,k)= CB_index_Y(i,j,k) * RA_index_Y_z(i,j,k);  
                
             CB_RB_Y_x(i,j,k)= CB_index_Y(i,j,k) *  RB_index_Y_x(i,j,k);  
             CB_RB_Y_z(i,j,k)= CB_index_Y(i,j,k) *  RB_index_Y_z(i,j,k);  
             
             CB_RA_sigma_s_Y_Y(i,j,k)= sigma_s_index_Y(i,j,k) * CB_index_Y(i,j,k) * RA_index_Y(i,j,k);  
             CB_RB_sigma_s_Y_Y(i,j,k)= sigma_s_index_Y(i,j,k) * CB_index_Y(i,j,k) * RB_index_Y(i,j,k);
             
         end
     end
 end
 
CB_RA_Z_y = zeros(Imax+1, Jmax+1, Kmax+1);
CB_RA_Z_x = zeros(Imax+1, Jmax+1, Kmax+1);
                
CB_RB_Z_y = zeros(Imax+1, Jmax+1, Kmax+1);
CB_RB_Z_x = zeros(Imax+1, Jmax+1, Kmax+1);
             
CB_RA_sigma_s_Z_Z = zeros(Imax+1, Jmax+1, Kmax+1);
CB_RB_sigma_s_Z_Z = zeros(Imax+1, Jmax+1, Kmax+1);

 for i=1:Imax+1
	 for j=1:Jmax+1
	     for k=1:Kmax+1
                
             CB_RA_Z_y(i,j,k)= CB_index_Z(i,j,k) * RA_index_Z_y(i,j,k);  
             CB_RA_Z_x(i,j,k)= CB_index_Z(i,j,k) * RA_index_Z_x(i,j,k);  
                
             CB_RB_Z_y(i,j,k)= CB_index_Z(i,j,k) *  RB_index_Z_y(i,j,k);  
             CB_RB_Z_x(i,j,k)= CB_index_Z(i,j,k) *  RB_index_Z_x(i,j,k);  
             
             CB_RA_sigma_s_Z_Z(i,j,k)= sigma_s_index_Z(i,j,k) * CB_index_Z(i,j,k) * RA_index_Z(i,j,k);  
             CB_RB_sigma_s_Z_Z(i,j,k)= sigma_s_index_Z(i,j,k) * CB_index_Z(i,j,k) * RB_index_Z(i,j,k);
             
         end
     end
 end 
 
v1_Ex = zeros(Imax+1, Jmax+1, Kmax+1);
v2_Ex = zeros(Imax+1, Jmax+1, Kmax+1); 
 
v1_Ey = zeros(Imax+1, Jmax+1, Kmax+1);
v2_Ey = zeros(Imax+1, Jmax+1, Kmax+1); 

v1_Ez = zeros(Imax+1, Jmax+1, Kmax+1);
v2_Ez = zeros(Imax+1, Jmax+1, Kmax+1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RAH_index_X_y=zeros(Imax+1, Jmax+1, Kmax+1);
RBH_index_X_y=zeros(Imax+1, Jmax+1, Kmax+1);
REH_index_X_y=zeros(Imax+1, Jmax+1, Kmax+1);
RFH_index_X_y=zeros(Imax+1, Jmax+1, Kmax+1);

RAH_index_X_z=zeros(Imax+1, Jmax+1, Kmax+1);
RBH_index_X_z=zeros(Imax+1, Jmax+1, Kmax+1);
REH_index_X_z=zeros(Imax+1, Jmax+1, Kmax+1);
RFH_index_X_z=zeros(Imax+1, Jmax+1, Kmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RAH_index_Y_x=zeros(Imax+1, Jmax+1, Kmax+1);
RBH_index_Y_x=zeros(Imax+1, Jmax+1, Kmax+1);
REH_index_Y_x=zeros(Imax+1, Jmax+1, Kmax+1);
RFH_index_Y_x=zeros(Imax+1, Jmax+1, Kmax+1);

RAH_index_Y_z=zeros(Imax+1, Jmax+1, Kmax+1);
RBH_index_Y_z=zeros(Imax+1, Jmax+1, Kmax+1);
REH_index_Y_z=zeros(Imax+1, Jmax+1, Kmax+1);
RFH_index_Y_z=zeros(Imax+1, Jmax+1, Kmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RAH_index_Z_x=zeros(Imax+1, Jmax+1, Kmax+1);
RBH_index_Z_x=zeros(Imax+1, Jmax+1, Kmax+1);
REH_index_Z_x=zeros(Imax+1, Jmax+1, Kmax+1);
RFH_index_Z_x=zeros(Imax+1, Jmax+1, Kmax+1);

RAH_index_Z_y=zeros(Imax+1, Jmax+1, Kmax+1);
RBH_index_Z_y=zeros(Imax+1, Jmax+1, Kmax+1);
REH_index_Z_y=zeros(Imax+1, Jmax+1, Kmax+1);
RFH_index_Z_y=zeros(Imax+1, Jmax+1, Kmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
			RAH_index_X_y(i,j,k) = RAh_y(j);
			RAH_index_X_z(i,j,k) = RAh_z(k);

			RBH_index_X_y(i,j,k) = RBh_y(j);
			RBH_index_X_z(i,j,k) = RBh_z(k);

			REH_index_X_y(i,j,k) = REh_y(j);
			REH_index_X_z(i,j,k) = REh_z(k);

			RFH_index_X_y(i,j,k) = RFh_y(j);
			RFH_index_X_z(i,j,k) = RFh_z(k);
        end
    end
end

for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
			RAH_index_Y_x(i,j,k) = RAh_x(i);
			RAH_index_Y_z(i,j,k) = RAh_z(k);

			RBH_index_Y_x(i,j,k) = RBh_x(i);
			RBH_index_Y_z(i,j,k) = RBh_z(k);

			REH_index_Y_x(i,j,k) = REh_x(i);
			REH_index_Y_z(i,j,k) = REh_z(k);

			RFH_index_Y_x(i,j,k) = RFh_x(i);
			RFH_index_Y_z(i,j,k) = RFh_z(k);
        end
    end
end

for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
			RAH_index_Z_x(i,j,k) = RAh_x(i);
			RAH_index_Z_y(i,j,k) = RAh_y(j);

			RBH_index_Z_x(i,j,k) = RBh_x(i);
			RBH_index_Z_y(i,j,k) = RBh_y(j);

			REH_index_Z_x(i,j,k) = REh_x(i);
			REH_index_Z_y(i,j,k) = REh_y(j);

			RFH_index_Z_x(i,j,k) = RFh_x(i);
			RFH_index_Z_y(i,j,k) = RFh_y(j);
        end
    end
end

v1_Hx = zeros(Imax+1, Jmax+1, Kmax+1);
v2_Hx = zeros(Imax+1, Jmax+1, Kmax+1);

v1_Hy = zeros(Imax+1, Jmax+1, Kmax+1);
v2_Hy = zeros(Imax+1, Jmax+1, Kmax+1);

v1_Hz = zeros(Imax+1, Jmax+1, Kmax+1);
v2_Hz = zeros(Imax+1, Jmax+1, Kmax+1);
%%  仿真计时
gpucompute_RI;
tic;
%%
%  BEGIN TIME STEP 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
for n = 1:Ntimesteps
    n
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%  UPDATE Ex 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    i=1:Imax;
	j=2:Jmax;
	k=2:Kmax;
            
	v1_Ex(i,j,k) = (Hz(i,j,k) - Hz(i,j-1,k)) / dy;
	v2_Ex(i,j,k) = (Hy(i,j,k) - Hy(i,j,k-1)) / dz;
              
	Ex_1(i,j,k) = Ex(i,j,k);
	Ex(i,j,k) = CA_Tx(i,j,k) .* Ex(i,j,k) + CB_Tx(i,j,k) .* (v1_Ex(i,j,k) - v2_Ex(i,j,k)) ...
              + CB_RA_X_y(i,j,k) .* v1_Ex(i,j,k) - CB_RA_X_z(i,j,k) .* v2_Ex(i,j,k)...
			  + CB_RB_X_y(i,j,k) .* GEZ_Sxy(i,j,k)...
		      - CB_RB_X_z(i,j,k) .* GEZ_Sxz(i,j,k) ...
			  - CB_RA_sigma_s_X_X(i,j,k) .* Ex_1(i,j,k) - CB_RB_sigma_s_X_X(i,j,k) .* GP_Sx(i,j,k);
		  
    % 辅助变量
	GEZ_Sxy(i,j,k) = RE_index_X_y(i,j,k) .* GEZ_Sxy(i,j,k) + RF_index_X_y(i,j,k) .* v1_Ex(i,j,k);
	GEZ_Sxz(i,j,k) = RE_index_X_z(i,j,k) .* GEZ_Sxz(i,j,k) + RF_index_X_z(i,j,k) .* v2_Ex(i,j,k);
	GP_Sx(i,j,k) = RE_index_X(i,j,k) .* GP_Sx(i,j,k) + RF_index_X(i,j,k) .* Ex_1(i,j,k);
      
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%  UPDATE Ey 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    j=1:Jmax;
	i=2:Imax;
    k=2:Kmax;

	v1_Ey(i,j,k) = (Hx(i,j,k) - Hx(i,j,k-1)) / dz;
	v2_Ey(i,j,k) = (Hz(i,j,k) - Hz(i-1,j,k)) / dx;

	Ey_1(i,j,k) = Ey(i,j,k);
	Ey(i,j,k) = CA_Ty(i,j,k) .* Ey(i,j,k) + CB_Ty(i,j,k) .* (v1_Ey(i,j,k) - v2_Ey(i,j,k))...
              + CB_RA_Y_z(i,j,k) .* v1_Ey(i,j,k) - CB_RA_Y_x(i,j,k) .* v2_Ey(i,j,k)...
			  + CB_RB_Y_z(i,j,k) .* GEZ_Syz(i,j,k)...
			  - CB_RB_Y_x(i,j,k) .* GEZ_Syx(i,j,k)...
			  - CB_RA_sigma_s_Y_Y(i,j,k) .* Ey_1(i,j,k) - CB_RB_sigma_s_Y_Y(i,j,k) .* GP_Sy(i,j,k);

	% 辅助变量
	GEZ_Syz(i,j,k) = RE_index_Y_z(i,j,k) .* GEZ_Syz(i,j,k) + RF_index_Y_z(i,j,k) .* v1_Ey(i,j,k);
	GEZ_Syx(i,j,k) = RE_index_Y_x(i,j,k) .* GEZ_Syx(i,j,k) + RF_index_Y_x(i,j,k) .* v2_Ey(i,j,k);
	GP_Sy(i,j,k) = RE_index_Y(i,j,k) .* GP_Sy(i,j,k) + RF_index_Y(i,j,k) .* Ey_1(i,j,k);
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%  UPDATE Ez
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    k=1:Kmax;
	i=2:Imax;
	j=2:Jmax;
            
	v1_Ez(i,j,k) = (Hy(i,j,k) - Hy(i-1,j,k)) / dx;
	v2_Ez(i,j,k) = (Hx(i,j,k) - Hx(i,j-1,k)) / dy;

	Ez_1(i,j,k) = Ez(i,j,k);
	Ez(i,j,k) = CA_Tz(i,j,k) .* Ez(i,j,k) + CB_Tz(i,j,k) .* (v1_Ez(i,j,k) - v2_Ez(i,j,k))...
			  + CB_RA_Z_x(i,j,k) .* v1_Ez(i,j,k) - CB_RA_Z_y(i,j,k) .* v2_Ez(i,j,k)...
			  + CB_RB_Z_x(i,j,k) .* GEZ_Szx(i,j,k)...
			  - CB_RB_Z_y(i,j,k) .* GEZ_Szy(i,j,k)...
			  - CB_RA_sigma_s_Z_Z(i,j,k) .* Ez_1(i,j,k) - CB_RB_sigma_s_Z_Z(i,j,k) .* GP_Sz(i,j,k);
    
	% 辅助变量
	GEZ_Szx(i,j,k) = RE_index_Z_x(i,j,k) .* GEZ_Szx(i,j,k) + RF_index_Z_x(i,j,k) .* v1_Ez(i,j,k);
	GEZ_Szy(i,j,k) = RE_index_Z_y(i,j,k) .* GEZ_Szy(i,j,k) + RF_index_Z_y(i,j,k) .* v2_Ez(i,j,k);
	GP_Sz(i,j,k) = RE_index_Z(i,j,k) .* GP_Sz(i,j,k) + RF_index_Z(i,j,k) .* Ez_1(i,j,k);

%%   源
%   source(n) = exp(-(((n + 1) * EM_DELTAT - t0) / (SPREAD))^2);
%   Ez(isource01+0,jsource01-0,ksource01)= Ez(isource01+0,jsource01-0,ksource01) + source(n);
    source(n) = 100 * exp(-(((n + 1 - 1) * EM_DELTAT - T0) / (SPREAD))^2);
%     source(n) = -2.0 * (((n + 1 - 1) * EM_DELTAT - T0) / (SPREAD)) * exp(-(((n + 1 - 1) * EM_DELTAT - T0) / (SPREAD))^2);
    %      Ez(IC,JC-10,KC) = source(n);
    Ez(IC,JC - 162,KC) = source(n);
% Probe Point   
    A_RI_1(n) = Ez(IC,JC - 12,KC);
    A_RI_2(n) = Ez(IC,JC + R + 12,KC);
      
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%  UPDATE Hx 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  i=2:Imax;
  j=1:Jmax;
  k=1:Kmax;
  % Hx Hxy Hxz
  v1_Hx(i,j,k) = (Ez(i,j+1,k) - Ez(i,j,k)) / dy;
  v2_Hx(i,j,k) = (Ey(i,j,k+1) - Ey(i,j,k)) / dz;

  Hx(i,j,k) = Hx(i,j,k) + CQ * (-v1_Hx(i,j,k) + v2_Hx(i,j,k))...
					    - CQ * (RAH_index_X_y(i,j,k) .* v1_Hx(i,j,k) + RBH_index_X_y(i,j,k) .* Hxy(i,j,k))...
					    + CQ * (RAH_index_X_z(i,j,k) .* v2_Hx(i,j,k) + RBH_index_X_z(i,j,k) .* Hxz(i,j,k));

  Hxy(i,j,k) = REH_index_X_y(i,j,k) .* Hxy(i,j,k) + RFH_index_X_y(i,j,k) .* v1_Hx(i,j,k);
  Hxz(i,j,k) = REH_index_X_z(i,j,k) .* Hxz(i,j,k) + RFH_index_X_z(i,j,k) .* v2_Hx(i,j,k);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%  UPDATE Hy 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  j=2:Jmax;
  i=1:Imax;
  k=1:Kmax;
  % Hy Hyx Hyz
  v1_Hy(i,j,k) = (Ex(i,j,k+1) - Ex(i,j,k)) / dz;
  v2_Hy(i,j,k) = (Ez(i+1,j,k) - Ez(i,j,k)) / dx;

  Hy(i,j,k) = Hy(i,j,k) + CQ * (-v1_Hy(i,j,k) + v2_Hy(i,j,k))...
					    - CQ * (RAH_index_Y_z(i,j,k) .* v1_Hy(i,j,k) + RBH_index_Y_z(i,j,k) .* Hyz(i,j,k))...
						+ CQ * (RAH_index_Y_x(i,j,k) .* v2_Hy(i,j,k) + RBH_index_Y_x(i,j,k) .* Hyx(i,j,k));

  Hyz(i,j,k) = REH_index_Y_z(i,j,k) .* Hyz(i,j,k) + RFH_index_Y_z(i,j,k) .* v1_Hy(i,j,k);
  Hyx(i,j,k) = REH_index_Y_x(i,j,k) .* Hyx(i,j,k) + RFH_index_Y_x(i,j,k) .* v2_Hy(i,j,k);
  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%  UPDATE Hz 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  k=2:Kmax;
  i=1:Imax;
  j=1:Jmax;
  % Hz Hzx Hzy
  v1_Hz(i,j,k) = (Ey(i+1,j,k) - Ey(i,j,k)) / dx;
  v2_Hz(i,j,k) = (Ex(i,j+1,k) - Ex(i,j,k)) / dy;

  Hz(i,j,k) = Hz(i,j,k) + CQ * (-v1_Hz(i,j,k) + v2_Hz(i,j,k))...
					    - CQ * (RAH_index_Z_x(i,j,k) .* v1_Hz(i,j,k) + RBH_index_Z_x(i,j,k) .* Hzx(i,j,k))...
						+ CQ * (RAH_index_Z_y(i,j,k) .* v2_Hz(i,j,k) + RBH_index_Z_y(i,j,k) .* Hzy(i,j,k));

  Hzx(i,j,k) = REH_index_Z_x(i,j,k) .* Hzx(i,j,k) + RFH_index_Z_x(i,j,k) .* v1_Hz(i,j,k);
  Hzy(i,j,k) = REH_index_Z_y(i,j,k) .* Hzy(i,j,k) + RFH_index_Z_y(i,j,k) .* v2_Hz(i,j,k);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%   VISUALIZATION 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

%{
if mod(n,N)==0; 
   timestep=int2str(n); 
   mesh(Ez(:,:,ksource01+5)); 
   axis([0 Imax 0 Jmax -5e-4 5e-4]); 
%    title(['Ez at time step =',timestep]) 
   pause(0.0001) 
end
%}


end
toc;