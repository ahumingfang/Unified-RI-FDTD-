clc;
clear; 
close all;  
 
fp = fopen('Head_voxel.dat');
data=fread(fp,'int8');
Brain_data=reshape(data,256,256,128);
Brain_index=(Brain_data(:,:,64));

% Grid for EM fields
Imax = 298;  
Jmax = 298; 
Kmax = 170; 

Ntimesteps =9347; 
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
dx = 1E-3;
dy = 1E-3; 
dz = 1E-3; 

% CFL Limit
CFL=0.5;
EM_DELTAT =  CFL * 1 / c0 / (((1 /dx/dx + 1 /dy/dy + 1 /dz/dz))^0.5); 
       dt =  CFL * 1 / c0 / (((1 /dx/dx + 1 /dy/dy + 1 /dz/dz))^0.5); 
       
SPREAD = 0.5*1.0618e-9;                 % pulse width  ***  5.0*1.e-15*ref_s (narrow)  15.e-16
T0 = 6*SPREAD;                % pulse center ***

% Computational Center            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isource01 = round(Imax/2);  
jsource01 = round(Jmax/2); 
ksource01 = round(Kmax/2); 

% PML thickness in each direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NPML = 11; 

%%%%%%%%%%%%%%%%%%%%   building index
Brain_p=zeros(Imax+1, Jmax+1, Kmax+1);

for i = 1:256
    for j = 1:256
        for k = 1:128
            Brain_p(i+NPML+5,j+NPML+5,k+NPML+10) =  Brain_data(i,j,k);
        end
    end
end

% PML Parameters
alpha = 4.0; 
alpha_aa = 1.0; 
 

max_sigma_x = 1.0 * (alpha + 1) / (150 * pi * dx);
max_sigma_y = 1.0 * (alpha + 1) / (150 * pi * dy);
max_sigma_z = 1.0 * (alpha + 1) / (150 * pi * dz); 

max_alpha_x = 0.005; 
max_alpha_y = 0.005; 
max_alpha_z = 0.005;

 
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

CA = 1.0; 
CB = (dt/(eps0)); 
 
CP = 1.0; 
CQ = (dt/(mu0)); 
  
                                                                 
% source location

IC = round((Imax + 1) / 2.0) + 1;
JC = round((Jmax + 1) / 2.0) + 1;
KC = round((Kmax + 1) / 2.0) + 1;

% Parameters of Debye model for dispaesive media

% (Debye_1)
% skin
e0_inf_Debye_1 = 29.85054779;
es_Debye_1 = 47.93014336;
tao_Debye_1 = 43.6257593e-12;
deltae_Debye_1 = es_Debye_1 - e0_inf_Debye_1;
sigma_Debye_1_S = 0.54073572;

CA_Debye_1 = (eps0 / dt - sigma_Debye_1_S / 2.0) / (eps0 / dt + sigma_Debye_1_S / 2.0);
CB_Debye_1 = 1.0 / (eps0 / dt + sigma_Debye_1_S / 2.0);

K_Debye_1 = e0_inf_Debye_1;
sigma_Debye_1 = deltae_Debye_1;
Alfa_Debye_1 = 1.0;

RA_Debye_1 = ((2.0 * tao_Debye_1 + dt * Alfa_Debye_1) / (2.0 * tao_Debye_1 * K_Debye_1 + dt * (Alfa_Debye_1 * K_Debye_1 + sigma_Debye_1))) - 1.0;
RB_Debye_1 = (2.0 * tao_Debye_1 * K_Debye_1) / (2.0 * tao_Debye_1 * K_Debye_1 + dt * (Alfa_Debye_1 * K_Debye_1 + sigma_Debye_1));
RC_Debye_1 = dt * (Alfa_Debye_1 - Alfa_Debye_1 * K_Debye_1 - sigma_Debye_1) / (tao_Debye_1 * K_Debye_1);
RD_Debye_1 = dt * (Alfa_Debye_1 * K_Debye_1 + sigma_Debye_1) / (tao_Debye_1 * K_Debye_1);
RE_Debye_1 = 1.0 - RD_Debye_1 * RB_Debye_1;
RF_Debye_1 = RC_Debye_1 - RA_Debye_1 * RD_Debye_1;

% (Debye_2)
% cerebrospinal fluid 
e0_inf_Debye_2 = 33.14797211;
es_Debye_2 = 70.39963913;
tao_Debye_2 = 18.1622616e-12;
deltae_Debye_2 = es_Debye_2 - e0_inf_Debye_2;
sigma_Debye_2_S = 2.14390564;

CA_Debye_2 = (eps0 / dt - sigma_Debye_2_S / 2.0) / (eps0 / dt + sigma_Debye_2_S / 2.0);
CB_Debye_2 = 1.0 / (eps0 / dt + sigma_Debye_2_S / 2.0);

K_Debye_2 = e0_inf_Debye_2;
sigma_Debye_2 = deltae_Debye_2;
Alfa_Debye_2 = 1.0;

RA_Debye_2 = ((2.0 * tao_Debye_2 + dt * Alfa_Debye_2) / (2.0 * tao_Debye_2 * K_Debye_2 + dt * (Alfa_Debye_2 * K_Debye_2 + sigma_Debye_2))) - 1.0;
RB_Debye_2 = (2.0 * tao_Debye_2 * K_Debye_2) / (2.0 * tao_Debye_2 * K_Debye_2 + dt * (Alfa_Debye_2 * K_Debye_2 + sigma_Debye_2));
RC_Debye_2 = dt * (Alfa_Debye_2 - Alfa_Debye_2 * K_Debye_2 - sigma_Debye_2) / (tao_Debye_2 * K_Debye_2);
RD_Debye_2 = dt * (Alfa_Debye_2 * K_Debye_2 + sigma_Debye_2) / (tao_Debye_2 * K_Debye_2);
RE_Debye_2 = 1.0 - RD_Debye_2 * RB_Debye_2;
RF_Debye_2 = RC_Debye_2 - RA_Debye_2 * RD_Debye_2;

% (Debye_3)
% spinal cord
e0_inf_Debye_3 = 20.98018456;
es_Debye_3 = 34.74881363;
tao_Debye_3 = 36.5079009e-12;
deltae_Debye_3 = es_Debye_3 - e0_inf_Debye_3;
sigma_Debye_3_S = 0.35986477;

CA_Debye_3 = (eps0 / dt - sigma_Debye_3_S / 2.0) / (eps0 / dt + sigma_Debye_3_S / 2.0);
CB_Debye_3 = 1.0 / (eps0 / dt + sigma_Debye_3_S / 2.0);

K_Debye_3 = e0_inf_Debye_3;
sigma_Debye_3 = deltae_Debye_3;
Alfa_Debye_3 = 1.0;

RA_Debye_3 = ((2.0 * tao_Debye_3 + dt * Alfa_Debye_3) / (2.0 * tao_Debye_3 * K_Debye_3 + dt * (Alfa_Debye_3 * K_Debye_3 + sigma_Debye_3))) - 1.0;
RB_Debye_3 = (2.0 * tao_Debye_3 * K_Debye_3) / (2.0 * tao_Debye_3 * K_Debye_3 + dt * (Alfa_Debye_3 * K_Debye_3 + sigma_Debye_3));
RC_Debye_3 = dt * (Alfa_Debye_3 - Alfa_Debye_3 * K_Debye_3 - sigma_Debye_3) / (tao_Debye_3 * K_Debye_3);
RD_Debye_3 = dt * (Alfa_Debye_3 * K_Debye_3 + sigma_Debye_3) / (tao_Debye_3 * K_Debye_3);
RE_Debye_3 = 1.0 - RD_Debye_3 * RB_Debye_3;
RF_Debye_3 = RC_Debye_3 - RA_Debye_3 * RD_Debye_3;

% (Debye_4)
% bone
e0_inf_Debye_4 = 7.36329556;
es_Debye_4 = 14.16905880;
tao_Debye_4 = 34.1064989e-12;
deltae_Debye_4 = es_Debye_4 - e0_inf_Debye_4;
sigma_Debye_4_S = 0.10401546;

CA_Debye_4 = (eps0 / dt - sigma_Debye_4_S / 2.0) / (eps0 / dt + sigma_Debye_4_S / 2.0);
CB_Debye_4 = 1.0 / (eps0 / dt + sigma_Debye_4_S / 2.0);

K_Debye_4 = e0_inf_Debye_4;
sigma_Debye_4 = deltae_Debye_4;
Alfa_Debye_4 = 1.0;

RA_Debye_4 = ((2.0 * tao_Debye_4 + dt * Alfa_Debye_4) / (2.0 * tao_Debye_4 * K_Debye_4 + dt * (Alfa_Debye_4 * K_Debye_4 + sigma_Debye_4))) - 1.0;
RB_Debye_4 = (2.0 * tao_Debye_4 * K_Debye_4) / (2.0 * tao_Debye_4 * K_Debye_4 + dt * (Alfa_Debye_4 * K_Debye_4 + sigma_Debye_4));
RC_Debye_4 = dt * (Alfa_Debye_4 - Alfa_Debye_4 * K_Debye_4 - sigma_Debye_4) / (tao_Debye_4 * K_Debye_4);
RD_Debye_4 = dt * (Alfa_Debye_4 * K_Debye_4 + sigma_Debye_4) / (tao_Debye_4 * K_Debye_4);
RE_Debye_4 = 1.0 - RD_Debye_4 * RB_Debye_4;
RF_Debye_4 = RC_Debye_4 - RA_Debye_4 * RD_Debye_4;

% (Debye_5)
% thyroid gland
e0_inf_Debye_5 = 27.73882103;
es_Debye_5 = 60.55047417;
tao_Debye_5 = 17.7273439e-12;
deltae_Debye_5 = es_Debye_5 - e0_inf_Debye_5;
sigma_Debye_5_S = 0.80898708;

CA_Debye_5 = (eps0 / dt - sigma_Debye_5_S / 2.0) / (eps0 / dt + sigma_Debye_5_S / 2.0);
CB_Debye_5 = 1.0 / (eps0 / dt + sigma_Debye_5_S / 2.0);

K_Debye_5 = e0_inf_Debye_5;
sigma_Debye_5 = deltae_Debye_5;
Alfa_Debye_5 = 1.0;

RA_Debye_5 = ((2.0 * tao_Debye_5 + dt * Alfa_Debye_5) / (2.0 * tao_Debye_5 * K_Debye_5 + dt * (Alfa_Debye_5 * K_Debye_5 + sigma_Debye_5))) - 1.0;
RB_Debye_5 = (2.0 * tao_Debye_5 * K_Debye_5) / (2.0 * tao_Debye_5 * K_Debye_5 + dt * (Alfa_Debye_5 * K_Debye_5 + sigma_Debye_5));
RC_Debye_5 = dt * (Alfa_Debye_5 - Alfa_Debye_5 * K_Debye_5 - sigma_Debye_5) / (tao_Debye_5 * K_Debye_5);
RD_Debye_5 = dt * (Alfa_Debye_5 * K_Debye_5 + sigma_Debye_5) / (tao_Debye_5 * K_Debye_5);
RE_Debye_5 = 1.0 - RD_Debye_5 * RB_Debye_5;
RF_Debye_5 = RC_Debye_5 - RA_Debye_5 * RD_Debye_5;

% (Debye_6)
% muscle
e0_inf_Debye_6 = 28.00134277;
es_Debye_6 = 56.93144607;
tao_Debye_6 = 18.6721420e-12;
deltae_Debye_6 = es_Debye_6 - e0_inf_Debye_6;
sigma_Debye_6_S = 0.74710494;

CA_Debye_6 = (eps0 / dt - sigma_Debye_6_S / 2.0) / (eps0 / dt + sigma_Debye_6_S / 2.0);
CB_Debye_6 = 1.0 / (eps0 / dt + sigma_Debye_6_S / 2.0);

K_Debye_6 = e0_inf_Debye_6;
sigma_Debye_6 = deltae_Debye_6;
Alfa_Debye_6 = 1.0;

RA_Debye_6 = ((2.0 * tao_Debye_6 + dt * Alfa_Debye_6) / (2.0 * tao_Debye_6 * K_Debye_6 + dt * (Alfa_Debye_6 * K_Debye_6 + sigma_Debye_6))) - 1.0;
RB_Debye_6 = (2.0 * tao_Debye_6 * K_Debye_6) / (2.0 * tao_Debye_6 * K_Debye_6 + dt * (Alfa_Debye_6 * K_Debye_6 + sigma_Debye_6));
RC_Debye_6 = dt * (Alfa_Debye_6 - Alfa_Debye_6 * K_Debye_6 - sigma_Debye_6) / (tao_Debye_6 * K_Debye_6);
RD_Debye_6 = dt * (Alfa_Debye_6 * K_Debye_6 + sigma_Debye_6) / (tao_Debye_6 * K_Debye_6);
RE_Debye_6 = 1.0 - RD_Debye_6 * RB_Debye_6;
RF_Debye_6 = RC_Debye_6 - RA_Debye_6 * RD_Debye_6;

% (Debye_7)
% Gland/Lacrimal glands (pineal gland )
e0_inf_Debye_7 = 27.70847893;
es_Debye_7 = 60.54281044;
tao_Debye_7 = 17.7068464e-12;
deltae_Debye_7 = es_Debye_7 - e0_inf_Debye_7;
sigma_Debye_7_S = 0.80907071;

CA_Debye_7 = (eps0 / dt - sigma_Debye_7_S / 2.0) / (eps0 / dt + sigma_Debye_7_S / 2.0);
CB_Debye_7 = 1.0 / (eps0 / dt + sigma_Debye_7_S / 2.0);

K_Debye_7 = e0_inf_Debye_7;
sigma_Debye_7 = deltae_Debye_7;
Alfa_Debye_7 = 1.0;

RA_Debye_7 = ((2.0 * tao_Debye_7 + dt * Alfa_Debye_7) / (2.0 * tao_Debye_7 * K_Debye_7 + dt * (Alfa_Debye_7 * K_Debye_7 + sigma_Debye_7))) - 1.0;
RB_Debye_7 = (2.0 * tao_Debye_7 * K_Debye_7) / (2.0 * tao_Debye_7 * K_Debye_7 + dt * (Alfa_Debye_7 * K_Debye_7 + sigma_Debye_7));
RC_Debye_7 = dt * (Alfa_Debye_7 - Alfa_Debye_7 * K_Debye_7 - sigma_Debye_7) / (tao_Debye_7 * K_Debye_7);
RD_Debye_7 = dt * (Alfa_Debye_7 * K_Debye_7 + sigma_Debye_7) / (tao_Debye_7 * K_Debye_7);
RE_Debye_7 = 1.0 - RD_Debye_7 * RB_Debye_7;
RF_Debye_7 = RC_Debye_7 - RA_Debye_7 * RD_Debye_7;

% (Debye_8)
% cerebellum 
e0_inf_Debye_8 = 35.19484711 ;
es_Debye_8 = 58.15469742;
tao_Debye_8 = 68.2758780e-12;
deltae_Debye_8 = es_Debye_8 - e0_inf_Debye_8;
sigma_Debye_8_S = 0.82605255 ;

CA_Debye_8 = (eps0 / dt - sigma_Debye_8_S / 2.0) / (eps0 / dt + sigma_Debye_8_S / 2.0);
CB_Debye_8 = 1.0 / (eps0 / dt + sigma_Debye_8_S / 2.0);

K_Debye_8 = e0_inf_Debye_8;
sigma_Debye_8 = deltae_Debye_8;
Alfa_Debye_8 = 1.0;

RA_Debye_8 = ((2.0 * tao_Debye_8 + dt * Alfa_Debye_8) / (2.0 * tao_Debye_8 * K_Debye_8 + dt * (Alfa_Debye_8 * K_Debye_8 + sigma_Debye_8))) - 1.0;
RB_Debye_8 = (2.0 * tao_Debye_8 * K_Debye_8) / (2.0 * tao_Debye_8 * K_Debye_8 + dt * (Alfa_Debye_8 * K_Debye_8 + sigma_Debye_8));
RC_Debye_8 = dt * (Alfa_Debye_8 - Alfa_Debye_8 * K_Debye_8 - sigma_Debye_8) / (tao_Debye_8 * K_Debye_8);
RD_Debye_8 = dt * (Alfa_Debye_8 * K_Debye_8 + sigma_Debye_8) / (tao_Debye_8 * K_Debye_8);
RE_Debye_8 = 1.0 - RD_Debye_8 * RB_Debye_8;
RF_Debye_8 = RC_Debye_8 - RA_Debye_8 * RD_Debye_8;

% (Debye_9)
% Tongue(Eye/Eyeball)
e0_inf_Debye_9 = 10.30803776;
es_Debye_9 = 67.71049214;
tao_Debye_9 = 8.27143302e-12;
deltae_Debye_9 = es_Debye_9 - e0_inf_Debye_9;
sigma_Debye_9_S = 1.44514167;

CA_Debye_9 = (eps0 / dt - sigma_Debye_9_S / 2.0) / (eps0 / dt + sigma_Debye_9_S / 2.0);
CB_Debye_9 = 1.0 / (eps0 / dt + sigma_Debye_9_S / 2.0);

K_Debye_9 = e0_inf_Debye_9;
sigma_Debye_9 = deltae_Debye_9;
Alfa_Debye_9 = 1.0;

RA_Debye_9 = ((2.0 * tao_Debye_9 + dt * Alfa_Debye_9) / (2.0 * tao_Debye_9 * K_Debye_9 + dt * (Alfa_Debye_9 * K_Debye_9 + sigma_Debye_9))) - 1.0;
RB_Debye_9 = (2.0 * tao_Debye_9 * K_Debye_9) / (2.0 * tao_Debye_9 * K_Debye_9 + dt * (Alfa_Debye_9 * K_Debye_9 + sigma_Debye_9));
RC_Debye_9 = dt * (Alfa_Debye_9 - Alfa_Debye_9 * K_Debye_9 - sigma_Debye_9) / (tao_Debye_9 * K_Debye_9);
RD_Debye_9 = dt * (Alfa_Debye_9 * K_Debye_9 + sigma_Debye_9) / (tao_Debye_9 * K_Debye_9);
RE_Debye_9 = 1.0 - RD_Debye_9 * RB_Debye_9;
RF_Debye_9 = RC_Debye_9 - RA_Debye_9 * RD_Debye_9;


% (Debye_11)  
% cartilage/Nasal septum
e0_inf_Debye_11 = 20.36329556;
es_Debye_11 = 38.583582;
tao_Debye_11 = 34.583582e-12;
deltae_Debye_11 = es_Debye_11 - e0_inf_Debye_11;
sigma_Debye_11_S =  0.49;

CA_Debye_11 = (eps0 / dt - sigma_Debye_11_S / 2.0) / (eps0 / dt + sigma_Debye_11_S / 2.0);
CB_Debye_11 = 1.0 / (eps0 / dt + sigma_Debye_11_S / 2.0);

K_Debye_11 = e0_inf_Debye_11;
sigma_Debye_11 = deltae_Debye_11;
Alfa_Debye_11 = 1.0;

RA_Debye_11 = ((2.0 * tao_Debye_11 + dt * Alfa_Debye_11) / (2.0 * tao_Debye_11 * K_Debye_11 + dt * (Alfa_Debye_11 * K_Debye_11 + sigma_Debye_11))) - 1.0;
RB_Debye_11 = (2.0 * tao_Debye_11 * K_Debye_11) / (2.0 * tao_Debye_11 * K_Debye_11 + dt * (Alfa_Debye_11 * K_Debye_11 + sigma_Debye_11));
RC_Debye_11 = dt * (Alfa_Debye_11 - Alfa_Debye_11 * K_Debye_11 - sigma_Debye_11) / (tao_Debye_11 * K_Debye_11);
RD_Debye_11 = dt * (Alfa_Debye_11 * K_Debye_11 + sigma_Debye_11) / (tao_Debye_11 * K_Debye_11);
RE_Debye_11 = 1.0 - RD_Debye_11 * RB_Debye_11;
RF_Debye_11 = RC_Debye_11 - RA_Debye_11 * RD_Debye_11;

% (Debye_12)
% white matter
e0_inf_Debye_12 = 24.37136650;
es_Debye_12 = 41.28079033;
tao_Debye_12 = 33.5850722e-12;
deltae_Debye_12 = es_Debye_12 - e0_inf_Debye_12;
sigma_Debye_12_S = 0.34836939;

CA_Debye_12 = (eps0 / dt - sigma_Debye_12_S / 2.0) / (eps0 / dt + sigma_Debye_12_S / 2.0);
CB_Debye_12 = 1.0 / (eps0 / dt + sigma_Debye_12_S / 2.0);

K_Debye_12 = e0_inf_Debye_12;
sigma_Debye_12 = deltae_Debye_12;
Alfa_Debye_12 = 1.0;

RA_Debye_12 = ((2.0 * tao_Debye_12 + dt * Alfa_Debye_12) / (2.0 * tao_Debye_12 * K_Debye_12 + dt * (Alfa_Debye_12 * K_Debye_12 + sigma_Debye_12))) - 1.0;
RB_Debye_12 = (2.0 * tao_Debye_12 * K_Debye_12) / (2.0 * tao_Debye_12 * K_Debye_12 + dt * (Alfa_Debye_12 * K_Debye_12 + sigma_Debye_12));
RC_Debye_12 = dt * (Alfa_Debye_12 - Alfa_Debye_12 * K_Debye_12 - sigma_Debye_12) / (tao_Debye_12 * K_Debye_12);
RD_Debye_12 = dt * (Alfa_Debye_12 * K_Debye_12 + sigma_Debye_12) / (tao_Debye_12 * K_Debye_12);
RE_Debye_12 = 1.0 - RD_Debye_12 * RB_Debye_12;
RF_Debye_12 = RC_Debye_12 - RA_Debye_12 * RD_Debye_12;

% (Debye_13)
% thalamus
e0_inf_Debye_13 = 33.05709076;
es_Debye_13 = 56.44398499;
tao_Debye_13 = 35.1972652e-12;
deltae_Debye_13 = es_Debye_13 - e0_inf_Debye_13;
sigma_Debye_13_S = 0.59515458 ;

CA_Debye_13 = (eps0 / dt - sigma_Debye_13_S / 2.0) / (eps0 / dt + sigma_Debye_13_S / 2.0);
CB_Debye_13 = 1.0 / (eps0 / dt + sigma_Debye_13_S / 2.0);

K_Debye_13 = e0_inf_Debye_13;
sigma_Debye_13 = deltae_Debye_13;
Alfa_Debye_13 = 1.0;

RA_Debye_13 = ((2.0 * tao_Debye_13 + dt * Alfa_Debye_13) / (2.0 * tao_Debye_13 * K_Debye_13 + dt * (Alfa_Debye_13 * K_Debye_13 + sigma_Debye_13))) - 1.0;
RB_Debye_13 = (2.0 * tao_Debye_13 * K_Debye_13) / (2.0 * tao_Debye_13 * K_Debye_13 + dt * (Alfa_Debye_13 * K_Debye_13 + sigma_Debye_13));
RC_Debye_13 = dt * (Alfa_Debye_13 - Alfa_Debye_13 * K_Debye_13 - sigma_Debye_13) / (tao_Debye_13 * K_Debye_13);
RD_Debye_13 = dt * (Alfa_Debye_13 * K_Debye_13 + sigma_Debye_13) / (tao_Debye_13 * K_Debye_13);
RE_Debye_13 = 1.0 - RD_Debye_13 * RB_Debye_13;
RF_Debye_13 = RC_Debye_13 - RA_Debye_13 * RD_Debye_13;

% (Debye_14) 
% Bone marrow （用spinal cord（脊髓）代替）
e0_inf_Debye_14 = 20.98018456 ;
es_Debye_14 = 34.74881363;
tao_Debye_14 = 36.5079009e-12;
deltae_Debye_14 = es_Debye_14 - e0_inf_Debye_14;
sigma_Debye_14_S = 0.35986477;

CA_Debye_14 = (eps0 / dt - sigma_Debye_14_S / 2.0) / (eps0 / dt + sigma_Debye_14_S / 2.0);
CB_Debye_14 = 1.0 / (eps0 / dt + sigma_Debye_14_S / 2.0);

K_Debye_14 = e0_inf_Debye_14;
sigma_Debye_14 = deltae_Debye_14;
Alfa_Debye_14 = 1.0;

RA_Debye_14 = ((2.0 * tao_Debye_14 + dt * Alfa_Debye_14) / (2.0 * tao_Debye_14 * K_Debye_14 + dt * (Alfa_Debye_14 * K_Debye_14 + sigma_Debye_14))) - 1.0;
RB_Debye_14 = (2.0 * tao_Debye_14 * K_Debye_14) / (2.0 * tao_Debye_14 * K_Debye_14 + dt * (Alfa_Debye_14 * K_Debye_14 + sigma_Debye_14));
RC_Debye_14 = dt * (Alfa_Debye_14 - Alfa_Debye_14 * K_Debye_14 - sigma_Debye_14) / (tao_Debye_14 * K_Debye_14);
RD_Debye_14 = dt * (Alfa_Debye_14 * K_Debye_14 + sigma_Debye_14) / (tao_Debye_14 * K_Debye_14);
RE_Debye_14 = 1.0 - RD_Debye_14 * RB_Debye_14;
RF_Debye_14 = RC_Debye_14 - RA_Debye_14 * RD_Debye_14;

% (Debye_15)
% pituitary gland 
e0_inf_Debye_15 = 27.70847893;
es_Debye_15 = 60.54281044;
tao_Debye_15 = 17.7068464e-12;
deltae_Debye_15 = es_Debye_15 - e0_inf_Debye_15;
sigma_Debye_15_S = 0.80907071;

CA_Debye_15 = (eps0 / dt - sigma_Debye_15_S / 2.0) / (eps0 / dt + sigma_Debye_15_S / 2.0);
CB_Debye_15 = 1.0 / (eps0 / dt + sigma_Debye_15_S / 2.0);

K_Debye_15 = e0_inf_Debye_15;
sigma_Debye_15 = deltae_Debye_15;
Alfa_Debye_15 = 1.0;

RA_Debye_15 = ((2.0 * tao_Debye_15 + dt * Alfa_Debye_15) / (2.0 * tao_Debye_15 * K_Debye_15 + dt * (Alfa_Debye_15 * K_Debye_15 + sigma_Debye_15))) - 1.0;
RB_Debye_15 = (2.0 * tao_Debye_15 * K_Debye_15) / (2.0 * tao_Debye_15 * K_Debye_15 + dt * (Alfa_Debye_15 * K_Debye_15 + sigma_Debye_15));
RC_Debye_15 = dt * (Alfa_Debye_15 - Alfa_Debye_15 * K_Debye_15 - sigma_Debye_15) / (tao_Debye_15 * K_Debye_15);
RD_Debye_15 = dt * (Alfa_Debye_15 * K_Debye_15 + sigma_Debye_15) / (tao_Debye_15 * K_Debye_15);
RE_Debye_15 = 1.0 - RD_Debye_15 * RB_Debye_15;
RF_Debye_15 = RC_Debye_15 - RA_Debye_15 * RD_Debye_15;

% (Debye_16)
% fat
e0_inf_Debye_16 = 3.99812841;
es_Debye_16 = 5.53071189;
tao_Debye_16 = 23.6329289e-12;
deltae_Debye_16 = es_Debye_16 - e0_inf_Debye_16;
sigma_Debye_16_S = 0.03710634;

CA_Debye_16 = (eps0 / dt - sigma_Debye_16_S / 2.0) / (eps0 / dt + sigma_Debye_16_S / 2.0);
CB_Debye_16 = 1.0 / (eps0 / dt + sigma_Debye_16_S / 2.0);

K_Debye_16 = e0_inf_Debye_16;
sigma_Debye_16 = deltae_Debye_16;
Alfa_Debye_16 = 1.0;

RA_Debye_16 = ((2.0 * tao_Debye_16 + dt * Alfa_Debye_16) / (2.0 * tao_Debye_16 * K_Debye_16 + dt * (Alfa_Debye_16 * K_Debye_16 + sigma_Debye_16))) - 1.0;
RB_Debye_16 = (2.0 * tao_Debye_16 * K_Debye_16) / (2.0 * tao_Debye_16 * K_Debye_16 + dt * (Alfa_Debye_16 * K_Debye_16 + sigma_Debye_16));
RC_Debye_16 = dt * (Alfa_Debye_16 - Alfa_Debye_16 * K_Debye_16 - sigma_Debye_16) / (tao_Debye_16 * K_Debye_16);
RD_Debye_16 = dt * (Alfa_Debye_16 * K_Debye_16 + sigma_Debye_16) / (tao_Debye_16 * K_Debye_16);
RE_Debye_16 = 1.0 - RD_Debye_16 * RB_Debye_16;
RF_Debye_16 = RC_Debye_16 - RA_Debye_16 * RD_Debye_16;

% (Debye_17)
% Blood pool(用资料中(cerebrospinal fluid(脑脊髓液)的参数设置))
e0_inf_Debye_17 = 33.14797211;
es_Debye_17 = 70.39963913;
tao_Debye_17 = 18.1622616e-12;
deltae_Debye_17 = es_Debye_17 - e0_inf_Debye_17;
sigma_Debye_17_S = 2.14390564;

CA_Debye_17 = (eps0 / dt - sigma_Debye_17_S / 2.0) / (eps0 / dt + sigma_Debye_17_S / 2.0);
CB_Debye_17 = 1.0 / (eps0 / dt + sigma_Debye_17_S / 2.0);

K_Debye_17 = e0_inf_Debye_17;
sigma_Debye_17 = deltae_Debye_17;
Alfa_Debye_17 = 1.0;

RA_Debye_17 = ((2.0 * tao_Debye_17 + dt * Alfa_Debye_17) / (2.0 * tao_Debye_17 * K_Debye_17 + dt * (Alfa_Debye_17 * K_Debye_17 + sigma_Debye_17))) - 1.0;
RB_Debye_17 = (2.0 * tao_Debye_17 * K_Debye_17) / (2.0 * tao_Debye_17 * K_Debye_17 + dt * (Alfa_Debye_17 * K_Debye_17 + sigma_Debye_17));
RC_Debye_17 = dt * (Alfa_Debye_17 - Alfa_Debye_17 * K_Debye_17 - sigma_Debye_17) / (tao_Debye_17 * K_Debye_17);
RD_Debye_17 = dt * (Alfa_Debye_17 * K_Debye_17 + sigma_Debye_17) / (tao_Debye_17 * K_Debye_17);
RE_Debye_17 = 1.0 - RD_Debye_17 * RB_Debye_17;
RF_Debye_17 = RC_Debye_17 - RA_Debye_17 * RD_Debye_17;

% (Debye_18)
% Cerebral falx (Dura) (用资料中(cerebral cortex(大脑皮层))的参数设置)
e0_inf_Debye_18 = 33.05709076 ;
es_Debye_18 = 56.44398499;
tao_Debye_18 = 35.1972652e-12;
deltae_Debye_18 = es_Debye_18 - e0_inf_Debye_18;
sigma_Debye_18_S = 0.59515458;

CA_Debye_18 = (eps0 / dt - sigma_Debye_18_S / 2.0) / (eps0 / dt + sigma_Debye_18_S / 2.0);
CB_Debye_18 = 1.0 / (eps0 / dt + sigma_Debye_18_S / 2.0);

K_Debye_18 = e0_inf_Debye_18;
sigma_Debye_18 = deltae_Debye_18;
Alfa_Debye_18 = 1.0;

RA_Debye_18 = ((2.0 * tao_Debye_18 + dt * Alfa_Debye_18) / (2.0 * tao_Debye_18 * K_Debye_18 + dt * (Alfa_Debye_18 * K_Debye_18 + sigma_Debye_18))) - 1.0;
RB_Debye_18 = (2.0 * tao_Debye_18 * K_Debye_18) / (2.0 * tao_Debye_18 * K_Debye_18 + dt * (Alfa_Debye_18 * K_Debye_18 + sigma_Debye_18));
RC_Debye_18 = dt * (Alfa_Debye_18 - Alfa_Debye_18 * K_Debye_18 - sigma_Debye_18) / (tao_Debye_18 * K_Debye_18);
RD_Debye_18 = dt * (Alfa_Debye_18 * K_Debye_18 + sigma_Debye_18) / (tao_Debye_18 * K_Debye_18);
RE_Debye_18 = 1.0 - RD_Debye_18 * RB_Debye_18;
RF_Debye_18 = RC_Debye_18 - RA_Debye_18 * RD_Debye_18;

% (Debye_19)
% Eye/Eye(cornea)
e0_inf_Debye_19 = 32.37210846;
es_Debye_19 = 57.84246635;
tao_Debye_19 = 28.7477611e-12;
deltae_Debye_19 = es_Debye_19 - e0_inf_Debye_19;
sigma_Debye_19_S = 1.06887341;

CA_Debye_19 = (eps0 / dt - sigma_Debye_19_S / 2.0) / (eps0 / dt + sigma_Debye_19_S / 2.0);
CB_Debye_19 = 1.0 / (eps0 / dt + sigma_Debye_19_S / 2.0);

K_Debye_19 = e0_inf_Debye_19;
sigma_Debye_19 = deltae_Debye_19;
Alfa_Debye_19 = 1.0;

RA_Debye_19 = ((2.0 * tao_Debye_19 + dt * Alfa_Debye_19) / (2.0 * tao_Debye_19 * K_Debye_19 + dt * (Alfa_Debye_19 * K_Debye_19 + sigma_Debye_19))) - 1.0;
RB_Debye_19 = (2.0 * tao_Debye_19 * K_Debye_19) / (2.0 * tao_Debye_19 * K_Debye_19 + dt * (Alfa_Debye_19 * K_Debye_19 + sigma_Debye_19));
RC_Debye_19 = dt * (Alfa_Debye_19 - Alfa_Debye_19 * K_Debye_19 - sigma_Debye_19) / (tao_Debye_19 * K_Debye_19);
RD_Debye_19 = dt * (Alfa_Debye_19 * K_Debye_19 + sigma_Debye_19) / (tao_Debye_19 * K_Debye_19);
RE_Debye_19 = 1.0 - RD_Debye_19 * RB_Debye_19;
RF_Debye_19 = RC_Debye_19 - RA_Debye_19 * RD_Debye_19;

% (Debye_20)
% Eye/Lens ((the crystalline lens) )
e0_inf_Debye_20 = 18.28276825;
es_Debye_20 = 37.82515335;
tao_Debye_20 = 21.0985136e-12;
deltae_Debye_20 = es_Debye_20 - e0_inf_Debye_20;
sigma_Debye_20_S = 0.34613180;

CA_Debye_20 = (eps0 / dt - sigma_Debye_20_S / 2.0) / (eps0 / dt + sigma_Debye_20_S / 2.0);
CB_Debye_20 = 1.0 / (eps0 / dt + sigma_Debye_20_S / 2.0);

K_Debye_20 = e0_inf_Debye_20;
sigma_Debye_20 = deltae_Debye_20;
Alfa_Debye_20 = 1.0;

RA_Debye_20 = ((2.0 * tao_Debye_20 + dt * Alfa_Debye_20) / (2.0 * tao_Debye_20 * K_Debye_20 + dt * (Alfa_Debye_20 * K_Debye_20 + sigma_Debye_20))) - 1.0;
RB_Debye_20 = (2.0 * tao_Debye_20 * K_Debye_20) / (2.0 * tao_Debye_20 * K_Debye_20 + dt * (Alfa_Debye_20 * K_Debye_20 + sigma_Debye_20));
RC_Debye_20 = dt * (Alfa_Debye_20 - Alfa_Debye_20 * K_Debye_20 - sigma_Debye_20) / (tao_Debye_20 * K_Debye_20);
RD_Debye_20 = dt * (Alfa_Debye_20 * K_Debye_20 + sigma_Debye_20) / (tao_Debye_20 * K_Debye_20);
RE_Debye_20 = 1.0 - RD_Debye_20 * RB_Debye_20;
RF_Debye_20 = RC_Debye_20 - RA_Debye_20 * RD_Debye_20;

%%%%%%%%%%%%%%
for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
		    CA_Tx(i,j,k) = CA;
		  	CB_Tx(i,j,k) = CB;
            
            if ( Brain_p(i,j,k) == 1)
                CA_Tx(i,j,k) = CA_Debye_1;
				CB_Tx(i,j,k) = CB_Debye_1;
                
            elseif ( Brain_p(i,j,k) == 2 || Brain_p(i,j,k) == 75 || Brain_p(i,j,k) == 92 || Brain_p(i,j,k) == 115 || Brain_p(i,j,k) == 122 || Brain_p(i,j,k) == 123)
				CA_Tx(i,j,k) = CA_Debye_2;
				CB_Tx(i,j,k) = CB_Debye_2;
                
             elseif ( Brain_p(i,j,k) == 3 || Brain_p(i,j,k) == 106 || Brain_p(i,j,k) == 111)
				CA_Tx(i,j,k) = CA_Debye_3;
				CB_Tx(i,j,k) = CB_Debye_3;
            
             elseif ( Brain_p(i,j,k) == 4 || Brain_p(i,j,k) == 5 || Brain_p(i,j,k) == 70 || Brain_p(i,j,k) == 71 || Brain_p(i,j,k) == 76 || Brain_p(i,j,k) == 81 || Brain_p(i,j,k) == 99 || Brain_p(i,j,k) == 125)
				CA_Tx(i,j,k) = CA_Debye_4;
				CB_Tx(i,j,k) = CB_Debye_4;  
                
             elseif ( Brain_p(i,j,k) == 72)
				CA_Tx(i,j,k) = CA_Debye_5;
				CB_Tx(i,j,k) = CB_Debye_5;  
              
             elseif ( Brain_p(i,j,k) == 9 || Brain_p(i,j,k) == 102)
				CA_Tx(i,j,k) = CA_Debye_6;
				CB_Tx(i,j,k) = CB_Debye_6;    
                
             elseif ( Brain_p(i,j,k) == 74)
				CA_Tx(i,j,k) = CA_Debye_7;
				CB_Tx(i,j,k) = CB_Debye_7;  
                
             elseif ( Brain_p(i,j,k) == 77 || Brain_p(i,j,k) == 85 || Brain_p(i,j,k) == 91)
				CA_Tx(i,j,k) = CA_Debye_8;
				CB_Tx(i,j,k) = CB_Debye_8;   
                
             elseif ( Brain_p(i,j,k) == 78 || Brain_p(i,j,k) == 110)
				CA_Tx(i,j,k) = CA_Debye_9;
				CB_Tx(i,j,k) = CB_Debye_9;     
                 
                
             elseif ( Brain_p(i,j,k) == 82 || Brain_p(i,j,k) == 30 || Brain_p(i,j,k) == 100)
				CA_Tx(i,j,k) = CA_Debye_11;
				CB_Tx(i,j,k) = CB_Debye_11;    
                
             elseif ( Brain_p(i,j,k) == 83)
				CA_Tx(i,j,k) = CA_Debye_12;
				CB_Tx(i,j,k) = CB_Debye_12;    
                
             elseif ( Brain_p(i,j,k) == 89 || Brain_p(i,j,k) == 95 || Brain_p(i,j,k) == 96 || Brain_p(i,j,k) == 101 || Brain_p(i,j,k) == 103 || Brain_p(i,j,k) == 105 || Brain_p(i,j,k) == 107 ...
                   || Brain_p(i,j,k) == 108 || Brain_p(i,j,k) == 109 || Brain_p(i,j,k) == 112 || Brain_p(i,j,k) == 114 || Brain_p(i,j,k) == 117 || Brain_p(i,j,k) == 118 || Brain_p(i,j,k) == 120 || Brain_p(i,j,k) == 124)
				CA_Tx(i,j,k) = CA_Debye_13;
				CB_Tx(i,j,k) = CB_Debye_13;     
                
             elseif ( Brain_p(i,j,k) == 26)
				CA_Tx(i,j,k) = CA_Debye_14;
				CB_Tx(i,j,k) = CB_Debye_14;      
                
             elseif ( Brain_p(i,j,k) == 97)
				CA_Tx(i,j,k) = CA_Debye_15;
				CB_Tx(i,j,k) = CB_Debye_15;        
             
             elseif ( Brain_p(i,j,k) == 22 || Brain_p(i,j,k) == 98 || Brain_p(i,j,k) == 116)
				CA_Tx(i,j,k) = CA_Debye_16;
				CB_Tx(i,j,k) = CB_Debye_16;       
                
             elseif ( Brain_p(i,j,k) == 23)
				CA_Tx(i,j,k) = CA_Debye_17;
				CB_Tx(i,j,k) = CB_Debye_17;      
                
             elseif ( Brain_p(i,j,k) == 113)
				CA_Tx(i,j,k) = CA_Debye_18;
				CB_Tx(i,j,k) = CB_Debye_18;     
                
             elseif ( Brain_p(i,j,k) == 119)
				CA_Tx(i,j,k) = CA_Debye_19;
				CB_Tx(i,j,k) = CB_Debye_19;      
                
             elseif ( Brain_p(i,j,k) == 121)
				CA_Tx(i,j,k) = CA_Debye_20;
				CB_Tx(i,j,k) = CB_Debye_20;           
            end
            
        end
    end
end

%%%%%%%%%%%%%%
for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
		    CA_Ty(i,j,k) = CA;
		  	CB_Ty(i,j,k) = CB;
            
            if ( Brain_p(i,j,k) == 1)
                CA_Ty(i,j,k) = CA_Debye_1;
				CB_Ty(i,j,k) = CB_Debye_1;
                
            elseif ( Brain_p(i,j,k) == 2 || Brain_p(i,j,k) == 75 || Brain_p(i,j,k) == 92 || Brain_p(i,j,k) == 115 || Brain_p(i,j,k) == 122 || Brain_p(i,j,k) == 123)
				CA_Ty(i,j,k) = CA_Debye_2;
				CB_Ty(i,j,k) = CB_Debye_2;
                
             elseif ( Brain_p(i,j,k) == 3 || Brain_p(i,j,k) == 106 || Brain_p(i,j,k) == 111)
				CA_Ty(i,j,k) = CA_Debye_3;
				CB_Ty(i,j,k) = CB_Debye_3;
            
             elseif ( Brain_p(i,j,k) == 4 || Brain_p(i,j,k) == 5 || Brain_p(i,j,k) == 70 || Brain_p(i,j,k) == 71 || Brain_p(i,j,k) == 76 || Brain_p(i,j,k) == 81 || Brain_p(i,j,k) == 99 || Brain_p(i,j,k) == 125)
				CA_Ty(i,j,k) = CA_Debye_4;
				CB_Ty(i,j,k) = CB_Debye_4;  
                
             elseif ( Brain_p(i,j,k) == 72)
				CA_Ty(i,j,k) = CA_Debye_5;
				CB_Ty(i,j,k) = CB_Debye_5;  
              
             elseif ( Brain_p(i,j,k) == 9 || Brain_p(i,j,k) == 102)
				CA_Ty(i,j,k) = CA_Debye_6;
				CB_Ty(i,j,k) = CB_Debye_6;    
                
             elseif ( Brain_p(i,j,k) == 74)
				CA_Ty(i,j,k) = CA_Debye_7;
				CB_Ty(i,j,k) = CB_Debye_7;  
                
             elseif ( Brain_p(i,j,k) == 77 || Brain_p(i,j,k) == 85 || Brain_p(i,j,k) == 91)
				CA_Ty(i,j,k) = CA_Debye_8;
				CB_Ty(i,j,k) = CB_Debye_8;   
                
             elseif ( Brain_p(i,j,k) == 78 || Brain_p(i,j,k) == 110)
				CA_Ty(i,j,k) = CA_Debye_9;
				CB_Ty(i,j,k) = CB_Debye_9;     
                                
             elseif ( Brain_p(i,j,k) == 82 || Brain_p(i,j,k) == 30 || Brain_p(i,j,k) == 100)
				CA_Ty(i,j,k) = CA_Debye_11;
				CB_Ty(i,j,k) = CB_Debye_11;    
                
             elseif ( Brain_p(i,j,k) == 83)
				CA_Ty(i,j,k) = CA_Debye_12;
				CB_Ty(i,j,k) = CB_Debye_12;    
                
             elseif ( Brain_p(i,j,k) == 89 || Brain_p(i,j,k) == 95 || Brain_p(i,j,k) == 96 || Brain_p(i,j,k) == 101 || Brain_p(i,j,k) == 103 || Brain_p(i,j,k) == 105 || Brain_p(i,j,k) == 107 ...
                   || Brain_p(i,j,k) == 108 || Brain_p(i,j,k) == 109 || Brain_p(i,j,k) == 112 || Brain_p(i,j,k) == 114 || Brain_p(i,j,k) == 117 || Brain_p(i,j,k) == 118 || Brain_p(i,j,k) == 120 || Brain_p(i,j,k) == 124)
				CA_Ty(i,j,k) = CA_Debye_13;
				CB_Ty(i,j,k) = CB_Debye_13;     
                
             elseif ( Brain_p(i,j,k) == 26)
				CA_Ty(i,j,k) = CA_Debye_14;
				CB_Ty(i,j,k) = CB_Debye_14;      
                
             elseif ( Brain_p(i,j,k) == 97)
				CA_Ty(i,j,k) = CA_Debye_15;
				CB_Ty(i,j,k) = CB_Debye_15;        
             
             elseif ( Brain_p(i,j,k) == 22 || Brain_p(i,j,k) == 98 || Brain_p(i,j,k) == 116)
				CA_Ty(i,j,k) = CA_Debye_16;
				CB_Ty(i,j,k) = CB_Debye_16;       
                
             elseif ( Brain_p(i,j,k) == 23)
				CA_Ty(i,j,k) = CA_Debye_17;
				CB_Ty(i,j,k) = CB_Debye_17;      
                
             elseif ( Brain_p(i,j,k) == 113)
				CA_Ty(i,j,k) = CA_Debye_18;
				CB_Ty(i,j,k) = CB_Debye_18;     
                
             elseif ( Brain_p(i,j,k) == 119)
				CA_Ty(i,j,k) = CA_Debye_19;
				CB_Ty(i,j,k) = CB_Debye_19;      
                
             elseif ( Brain_p(i,j,k) == 121)
				CA_Ty(i,j,k) = CA_Debye_20;
				CB_Ty(i,j,k) = CB_Debye_20;           
            end
            
        end
    end
end

%%%%%%%%%%%%%%
for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
		    CA_Tz(i,j,k) = CA;
		  	CB_Tz(i,j,k) = CB;
            
            if ( Brain_p(i,j,k) == 1)
                CA_Tz(i,j,k) = CA_Debye_1;
				CB_Tz(i,j,k) = CB_Debye_1;
                
            elseif ( Brain_p(i,j,k) == 2 || Brain_p(i,j,k) == 75 || Brain_p(i,j,k) == 92 || Brain_p(i,j,k) == 115 || Brain_p(i,j,k) == 122 || Brain_p(i,j,k) == 123)
				CA_Tz(i,j,k) = CA_Debye_2;
				CB_Tz(i,j,k) = CB_Debye_2;
                
             elseif ( Brain_p(i,j,k) == 3 || Brain_p(i,j,k) == 106 || Brain_p(i,j,k) == 111)
				CA_Tz(i,j,k) = CA_Debye_3;
				CB_Tz(i,j,k) = CB_Debye_3;
            
             elseif ( Brain_p(i,j,k) == 4 || Brain_p(i,j,k) == 5 || Brain_p(i,j,k) == 70 || Brain_p(i,j,k) == 71 || Brain_p(i,j,k) == 76 || Brain_p(i,j,k) == 81 || Brain_p(i,j,k) == 99 || Brain_p(i,j,k) == 125)
				CA_Tz(i,j,k) = CA_Debye_4;
				CB_Tz(i,j,k) = CB_Debye_4;  
                
             elseif ( Brain_p(i,j,k) == 72)
				CA_Tz(i,j,k) = CA_Debye_5;
				CB_Tz(i,j,k) = CB_Debye_5;  
              
             elseif ( Brain_p(i,j,k) == 9 || Brain_p(i,j,k) == 102)
				CA_Tz(i,j,k) = CA_Debye_6;
				CB_Tz(i,j,k) = CB_Debye_6;    
                
             elseif ( Brain_p(i,j,k) == 74)
				CA_Tz(i,j,k) = CA_Debye_7;
				CB_Tz(i,j,k) = CB_Debye_7;  
                
             elseif ( Brain_p(i,j,k) == 77 || Brain_p(i,j,k) == 85 || Brain_p(i,j,k) == 91)
				CA_Tz(i,j,k) = CA_Debye_8;
				CB_Tz(i,j,k) = CB_Debye_8;   
                
             elseif ( Brain_p(i,j,k) == 78 || Brain_p(i,j,k) == 110)
				CA_Tz(i,j,k) = CA_Debye_9;
				CB_Tz(i,j,k) = CB_Debye_9;     
                                
             elseif ( Brain_p(i,j,k) == 82 || Brain_p(i,j,k) == 30 || Brain_p(i,j,k) == 100)
				CA_Tz(i,j,k) = CA_Debye_11;
				CB_Tz(i,j,k) = CB_Debye_11;    
                
             elseif ( Brain_p(i,j,k) == 83)
				CA_Tz(i,j,k) = CA_Debye_12;
				CB_Tz(i,j,k) = CB_Debye_12;    
                
             elseif ( Brain_p(i,j,k) == 89 || Brain_p(i,j,k) == 95 || Brain_p(i,j,k) == 96 || Brain_p(i,j,k) == 101 || Brain_p(i,j,k) == 103 || Brain_p(i,j,k) == 105 || Brain_p(i,j,k) == 107 ...
                   || Brain_p(i,j,k) == 108 || Brain_p(i,j,k) == 109 || Brain_p(i,j,k) == 112 || Brain_p(i,j,k) == 114 || Brain_p(i,j,k) == 117 || Brain_p(i,j,k) == 118 || Brain_p(i,j,k) == 120 || Brain_p(i,j,k) == 124)
				CA_Tz(i,j,k) = CA_Debye_13;
				CB_Tz(i,j,k) = CB_Debye_13;     
                
             elseif ( Brain_p(i,j,k) == 26)
				CA_Tz(i,j,k) = CA_Debye_14;
				CB_Tz(i,j,k) = CB_Debye_14;      
                
             elseif ( Brain_p(i,j,k) == 97)
				CA_Tz(i,j,k) = CA_Debye_15;
				CB_Tz(i,j,k) = CB_Debye_15;        
             
             elseif ( Brain_p(i,j,k) == 22 || Brain_p(i,j,k) == 98 || Brain_p(i,j,k) == 116)
				CA_Tz(i,j,k) = CA_Debye_16;
				CB_Tz(i,j,k) = CB_Debye_16;       
                
             elseif ( Brain_p(i,j,k) == 23)
				CA_Tz(i,j,k) = CA_Debye_17;
				CB_Tz(i,j,k) = CB_Debye_17;      
                
             elseif ( Brain_p(i,j,k) == 113)
				CA_Tz(i,j,k) = CA_Debye_18;
				CB_Tz(i,j,k) = CB_Debye_18;     
                
             elseif ( Brain_p(i,j,k) == 119)
				CA_Tz(i,j,k) = CA_Debye_19;
				CB_Tz(i,j,k) = CB_Debye_19;      
                
             elseif ( Brain_p(i,j,k) == 121)
				CA_Tz(i,j,k) = CA_Debye_20;
				CB_Tz(i,j,k) = CB_Debye_20;           
            end
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
index_rou=ones(Imax+1, Jmax+1, Kmax+1);          % Density rou (kg/(m^3))

for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
		    index_rou(i,j,k) = 1;
    
            if ( Brain_p(i,j,k) == 1)
                index_rou(i,j,k) = 1109;
                
            elseif ( Brain_p(i,j,k) == 2 || Brain_p(i,j,k) == 75 || Brain_p(i,j,k) == 92 || Brain_p(i,j,k) == 115 || Brain_p(i,j,k) == 122 || Brain_p(i,j,k) == 123)
				 index_rou(i,j,k) = 1007;
                
             elseif ( Brain_p(i,j,k) == 3)
				 index_rou(i,j,k) = 1850;
            
             elseif ( Brain_p(i,j,k) == 106 || Brain_p(i,j,k) == 111)
				 index_rou(i,j,k) = 1075;    
                 
             elseif ( Brain_p(i,j,k) == 4 || Brain_p(i,j,k) == 5 || Brain_p(i,j,k) == 70 || Brain_p(i,j,k) == 71 || Brain_p(i,j,k) == 76 || Brain_p(i,j,k) == 81 || Brain_p(i,j,k) == 99 || Brain_p(i,j,k) == 125)
				 index_rou(i,j,k) = 1908;
                
             elseif ( Brain_p(i,j,k) == 72)
				 index_rou(i,j,k) = 1048;
              
             elseif ( Brain_p(i,j,k) == 9 || Brain_p(i,j,k) == 102)
				 index_rou(i,j,k) = 1090;
                
             elseif ( Brain_p(i,j,k) == 74)
				 index_rou(i,j,k) = 1028;
                
             elseif ( Brain_p(i,j,k) == 77 || Brain_p(i,j,k) == 85 || Brain_p(i,j,k) == 91)
				 index_rou(i,j,k) = 1045;
                
             elseif ( Brain_p(i,j,k) == 78)
				 index_rou(i,j,k) = 1090;
                 
             elseif (Brain_p(i,j,k) == 110)
				 index_rou(i,j,k) = 1005;    
                
             elseif ( Brain_p(i,j,k) == 15 || Brain_p(i,j,k) == 84 || Brain_p(i,j,k) == 104)
				 index_rou(i,j,k) = 1;
                
             elseif ( Brain_p(i,j,k) == 82 || Brain_p(i,j,k) == 30 || Brain_p(i,j,k) == 100)
				 index_rou(i,j,k) = 1110;
                
             elseif ( Brain_p(i,j,k) == 83)
				 index_rou(i,j,k) = 1041;
                
             elseif ( Brain_p(i,j,k) == 89 || Brain_p(i,j,k) == 95 || Brain_p(i,j,k) == 96 || Brain_p(i,j,k) == 101 || Brain_p(i,j,k) == 103 || Brain_p(i,j,k) == 105 || Brain_p(i,j,k) == 107 ...
                   || Brain_p(i,j,k) == 108 || Brain_p(i,j,k) == 109 || Brain_p(i,j,k) == 112 || Brain_p(i,j,k) == 114 || Brain_p(i,j,k) == 117 || Brain_p(i,j,k) == 118 || Brain_p(i,j,k) == 120 || Brain_p(i,j,k) == 124)
				 index_rou(i,j,k) = 1045;
                
             elseif ( Brain_p(i,j,k) == 26)
				 index_rou(i,j,k) = 1029;
                
             elseif ( Brain_p(i,j,k) == 97)
				 index_rou(i,j,k) = 1050;   
             
             elseif ( Brain_p(i,j,k) == 22 || Brain_p(i,j,k) == 98 || Brain_p(i,j,k) == 116)
				 index_rou(i,j,k) = 911;
                
             elseif ( Brain_p(i,j,k) == 23)
				 index_rou(i,j,k) = 1050;
                
             elseif ( Brain_p(i,j,k) == 113)
				 index_rou(i,j,k) = 1174;
                
             elseif ( Brain_p(i,j,k) == 119)
				 index_rou(i,j,k) = 1005;
                
             elseif ( Brain_p(i,j,k) == 121)
				 index_rou(i,j,k) = 1076;
            end
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
index_sig=zeros(Imax+1, Jmax+1, Kmax+1);          % Density rou (kg/(m^3))

for i = 1:Imax + 1
	for j = 1:Jmax + 1
	    for k = 1:Kmax + 1
		    index_sig(i,j,k) = 0;
    
            if ( Brain_p(i,j,k) == 1)
                index_sig(i,j,k) = 0.54073572;
                
            elseif ( Brain_p(i,j,k) == 2 || Brain_p(i,j,k) == 75 || Brain_p(i,j,k) == 92 || Brain_p(i,j,k) == 115 || Brain_p(i,j,k) == 122 || Brain_p(i,j,k) == 123)
				 index_sig(i,j,k) = 2.14390564;
                
             elseif ( Brain_p(i,j,k) == 3)
				 index_sig(i,j,k) = 0.35986477;
            
             elseif ( Brain_p(i,j,k) == 106 || Brain_p(i,j,k) == 111)
				 index_sig(i,j,k) = 0.35986477;    
                 
             elseif ( Brain_p(i,j,k) == 4 || Brain_p(i,j,k) == 5 || Brain_p(i,j,k) == 70 || Brain_p(i,j,k) == 71 || Brain_p(i,j,k) == 76 || Brain_p(i,j,k) == 81 || Brain_p(i,j,k) == 99 || Brain_p(i,j,k) == 125)
				 index_sig(i,j,k) = 0.10401546;
                
             elseif ( Brain_p(i,j,k) == 72)
				 index_sig(i,j,k) = 0.80898708;
              
             elseif ( Brain_p(i,j,k) == 9 || Brain_p(i,j,k) == 102)
				 index_sig(i,j,k) = 0.74710494;
                
             elseif ( Brain_p(i,j,k) == 74)
				 index_sig(i,j,k) = 0.80907071;
                
             elseif ( Brain_p(i,j,k) == 77 || Brain_p(i,j,k) == 85 || Brain_p(i,j,k) == 91)
				 index_sig(i,j,k) = 0.82605255;
                
             elseif ( Brain_p(i,j,k) == 78)
				 index_sig(i,j,k) = 1.44514167;
                 
             elseif (Brain_p(i,j,k) == 110)
				 index_sig(i,j,k) = 1.44514167; 
                
             elseif ( Brain_p(i,j,k) == 15 || Brain_p(i,j,k) == 84 || Brain_p(i,j,k) == 104)
				 index_sig(i,j,k) = 0.00000000;
                
             elseif ( Brain_p(i,j,k) == 82 || Brain_p(i,j,k) == 30 || Brain_p(i,j,k) == 100)
				 index_sig(i,j,k) =  0.4900000;
                
             elseif ( Brain_p(i,j,k) == 83)
				 index_sig(i,j,k) = 0.34836939;
                
             elseif ( Brain_p(i,j,k) == 89 || Brain_p(i,j,k) == 95 || Brain_p(i,j,k) == 96 || Brain_p(i,j,k) == 101 || Brain_p(i,j,k) == 103 || Brain_p(i,j,k) == 105 || Brain_p(i,j,k) == 107 ...
                   || Brain_p(i,j,k) == 108 || Brain_p(i,j,k) == 109 || Brain_p(i,j,k) == 112 || Brain_p(i,j,k) == 114 || Brain_p(i,j,k) == 117 || Brain_p(i,j,k) == 118 || Brain_p(i,j,k) == 120 || Brain_p(i,j,k) == 124)
				 index_sig(i,j,k) = 0.59515458 ;
                
             elseif ( Brain_p(i,j,k) == 26)
				 index_sig(i,j,k) = 0.35986477;
                
             elseif ( Brain_p(i,j,k) == 97)
				 index_sig(i,j,k) = 0.80907071;
             
             elseif ( Brain_p(i,j,k) == 22 || Brain_p(i,j,k) == 98 || Brain_p(i,j,k) == 116)
				 index_sig(i,j,k) = 0.03710634;
                
             elseif ( Brain_p(i,j,k) == 23)
				 index_sig(i,j,k) = 2.14390564;
                
             elseif ( Brain_p(i,j,k) == 113)
				 index_sig(i,j,k) = 0.59515458;
                
             elseif ( Brain_p(i,j,k) == 119)
				 index_sig(i,j,k) = 1.06887341;
                
             elseif ( Brain_p(i,j,k) == 121)
				 index_sig(i,j,k) = 0.34613180;
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
            
            if ( Brain_p(i,j,k) == 1)
                
				CB_index_X(i,j,k) = CB_Debye_1;

				RA_index_X_y(i,j,k) = RA_Debye_1;
				RA_index_X_z(i,j,k) = RA_Debye_1;

				RB_index_X_y(i,j,k) = RB_Debye_1;
				RB_index_X_z(i,j,k) = RB_Debye_1;

				RE_index_X_y(i,j,k) = RE_Debye_1;
				RE_index_X_z(i,j,k) = RE_Debye_1;

				RF_index_X_y(i,j,k) = RF_Debye_1;
				RF_index_X_z(i,j,k) = RF_Debye_1;

				RA_index_X(i,j,k) = RA_Debye_1;
				RB_index_X(i,j,k) = RB_Debye_1;
				RE_index_X(i,j,k) = RE_Debye_1;
				RF_index_X(i,j,k) = RF_Debye_1;
				sigma_s_index_X(i,j,k) = sigma_Debye_1_S;
                
            elseif ( Brain_p(i,j,k) == 2 || Brain_p(i,j,k) == 75 || Brain_p(i,j,k) == 92 || Brain_p(i,j,k) == 115 || Brain_p(i,j,k) == 122 || Brain_p(i,j,k) == 123)
               
                CB_index_X(i,j,k) = CB_Debye_2;

				RA_index_X_y(i,j,k) = RA_Debye_2;
				RA_index_X_z(i,j,k) = RA_Debye_2;

				RB_index_X_y(i,j,k) = RB_Debye_2;
				RB_index_X_z(i,j,k) = RB_Debye_2;

				RE_index_X_y(i,j,k) = RE_Debye_2;
				RE_index_X_z(i,j,k) = RE_Debye_2;

				RF_index_X_y(i,j,k) = RF_Debye_2;
				RF_index_X_z(i,j,k) = RF_Debye_2;

				RA_index_X(i,j,k) = RA_Debye_2;
				RB_index_X(i,j,k) = RB_Debye_2;
				RE_index_X(i,j,k) = RE_Debye_2;
				RF_index_X(i,j,k) = RF_Debye_2;
				sigma_s_index_X(i,j,k) = sigma_Debye_2_S;
                
             elseif ( Brain_p(i,j,k) == 3 || Brain_p(i,j,k) == 106 || Brain_p(i,j,k) == 111)
				
                CB_index_X(i,j,k) = CB_Debye_3;

				RA_index_X_y(i,j,k) = RA_Debye_3;
				RA_index_X_z(i,j,k) = RA_Debye_3;

				RB_index_X_y(i,j,k) = RB_Debye_3;
				RB_index_X_z(i,j,k) = RB_Debye_3;

				RE_index_X_y(i,j,k) = RE_Debye_3;
				RE_index_X_z(i,j,k) = RE_Debye_3;

				RF_index_X_y(i,j,k) = RF_Debye_3;
				RF_index_X_z(i,j,k) = RF_Debye_3;

				RA_index_X(i,j,k) = RA_Debye_3;
				RB_index_X(i,j,k) = RB_Debye_3;
				RE_index_X(i,j,k) = RE_Debye_3;
				RF_index_X(i,j,k) = RF_Debye_3;
				sigma_s_index_X(i,j,k) = sigma_Debye_3_S;
                
             elseif ( Brain_p(i,j,k) == 4 || Brain_p(i,j,k) == 5 || Brain_p(i,j,k) == 70 || Brain_p(i,j,k) == 71 || Brain_p(i,j,k) == 76 || Brain_p(i,j,k) == 81 || Brain_p(i,j,k) == 99 || Brain_p(i,j,k) == 125)
			
                CB_index_X(i,j,k) = CB_Debye_4;

				RA_index_X_y(i,j,k) = RA_Debye_4;
				RA_index_X_z(i,j,k) = RA_Debye_4;

				RB_index_X_y(i,j,k) = RB_Debye_4;
				RB_index_X_z(i,j,k) = RB_Debye_4;

				RE_index_X_y(i,j,k) = RE_Debye_4;
				RE_index_X_z(i,j,k) = RE_Debye_4;

				RF_index_X_y(i,j,k) = RF_Debye_4;
				RF_index_X_z(i,j,k) = RF_Debye_4;

				RA_index_X(i,j,k) = RA_Debye_4;
				RB_index_X(i,j,k) = RB_Debye_4;
				RE_index_X(i,j,k) = RE_Debye_4;
				RF_index_X(i,j,k) = RF_Debye_4;
				sigma_s_index_X(i,j,k) = sigma_Debye_4_S;
                
             elseif ( Brain_p(i,j,k) == 72)
				
                CB_index_X(i,j,k) = CB_Debye_5;

				RA_index_X_y(i,j,k) = RA_Debye_5;
				RA_index_X_z(i,j,k) = RA_Debye_5;

				RB_index_X_y(i,j,k) = RB_Debye_5;
				RB_index_X_z(i,j,k) = RB_Debye_5;

				RE_index_X_y(i,j,k) = RE_Debye_5;
				RE_index_X_z(i,j,k) = RE_Debye_5;

				RF_index_X_y(i,j,k) = RF_Debye_5;
				RF_index_X_z(i,j,k) = RF_Debye_5;

				RA_index_X(i,j,k) = RA_Debye_5;
				RB_index_X(i,j,k) = RB_Debye_5;
				RE_index_X(i,j,k) = RE_Debye_5;
				RF_index_X(i,j,k) = RF_Debye_5;
				sigma_s_index_X(i,j,k) = sigma_Debye_5_S;
                  
             elseif ( Brain_p(i,j,k) == 9 || Brain_p(i,j,k) == 102)
				
                CB_index_X(i,j,k) = CB_Debye_6;

				RA_index_X_y(i,j,k) = RA_Debye_6;
				RA_index_X_z(i,j,k) = RA_Debye_6;

				RB_index_X_y(i,j,k) = RB_Debye_6;
				RB_index_X_z(i,j,k) = RB_Debye_6;

				RE_index_X_y(i,j,k) = RE_Debye_6;
				RE_index_X_z(i,j,k) = RE_Debye_6;

				RF_index_X_y(i,j,k) = RF_Debye_6;
				RF_index_X_z(i,j,k) = RF_Debye_6;

				RA_index_X(i,j,k) = RA_Debye_6;
				RB_index_X(i,j,k) = RB_Debye_6;
				RE_index_X(i,j,k) = RE_Debye_6;
				RF_index_X(i,j,k) = RF_Debye_6;
				sigma_s_index_X(i,j,k) = sigma_Debye_6_S;
                
             elseif ( Brain_p(i,j,k) == 74)
				
                CB_index_X(i,j,k) = CB_Debye_7;

				RA_index_X_y(i,j,k) = RA_Debye_7;
				RA_index_X_z(i,j,k) = RA_Debye_7;

				RB_index_X_y(i,j,k) = RB_Debye_7;
				RB_index_X_z(i,j,k) = RB_Debye_7;

				RE_index_X_y(i,j,k) = RE_Debye_7;
				RE_index_X_z(i,j,k) = RE_Debye_7;

				RF_index_X_y(i,j,k) = RF_Debye_7;
				RF_index_X_z(i,j,k) = RF_Debye_7;

				RA_index_X(i,j,k) = RA_Debye_7;
				RB_index_X(i,j,k) = RB_Debye_7;
				RE_index_X(i,j,k) = RE_Debye_7;
				RF_index_X(i,j,k) = RF_Debye_7;
				sigma_s_index_X(i,j,k) = sigma_Debye_7_S;
                
             elseif ( Brain_p(i,j,k) == 77 || Brain_p(i,j,k) == 85 || Brain_p(i,j,k) == 91)
				
                CB_index_X(i,j,k) = CB_Debye_8;

				RA_index_X_y(i,j,k) = RA_Debye_8;
				RA_index_X_z(i,j,k) = RA_Debye_8;

				RB_index_X_y(i,j,k) = RB_Debye_8;
				RB_index_X_z(i,j,k) = RB_Debye_8;

				RE_index_X_y(i,j,k) = RE_Debye_8;
				RE_index_X_z(i,j,k) = RE_Debye_8;

				RF_index_X_y(i,j,k) = RF_Debye_8;
				RF_index_X_z(i,j,k) = RF_Debye_8;

				RA_index_X(i,j,k) = RA_Debye_8;
				RB_index_X(i,j,k) = RB_Debye_8;
				RE_index_X(i,j,k) = RE_Debye_8;
				RF_index_X(i,j,k) = RF_Debye_8;
				sigma_s_index_X(i,j,k) = sigma_Debye_8_S;  
                
             elseif ( Brain_p(i,j,k) == 78 || Brain_p(i,j,k) == 110)
				
                CB_index_X(i,j,k) = CB_Debye_9;

				RA_index_X_y(i,j,k) = RA_Debye_9;
				RA_index_X_z(i,j,k) = RA_Debye_9;

				RB_index_X_y(i,j,k) = RB_Debye_9;
				RB_index_X_z(i,j,k) = RB_Debye_9;

				RE_index_X_y(i,j,k) = RE_Debye_9;
				RE_index_X_z(i,j,k) = RE_Debye_9;

				RF_index_X_y(i,j,k) = RF_Debye_9;
				RF_index_X_z(i,j,k) = RF_Debye_9;

				RA_index_X(i,j,k) = RA_Debye_9;
				RB_index_X(i,j,k) = RB_Debye_9;
				RE_index_X(i,j,k) = RE_Debye_9;
				RF_index_X(i,j,k) = RF_Debye_9;
				sigma_s_index_X(i,j,k) = sigma_Debye_9_S;
                                
             elseif ( Brain_p(i,j,k) == 82 || Brain_p(i,j,k) == 30 || Brain_p(i,j,k) == 100)
				
                CB_index_X(i,j,k) = CB_Debye_11;

				RA_index_X_y(i,j,k) = RA_Debye_11;
				RA_index_X_z(i,j,k) = RA_Debye_11;

				RB_index_X_y(i,j,k) = RB_Debye_11;
				RB_index_X_z(i,j,k) = RB_Debye_11;

				RE_index_X_y(i,j,k) = RE_Debye_11;
				RE_index_X_z(i,j,k) = RE_Debye_11;

				RF_index_X_y(i,j,k) = RF_Debye_11;
				RF_index_X_z(i,j,k) = RF_Debye_11;

				RA_index_X(i,j,k) = RA_Debye_11;
				RB_index_X(i,j,k) = RB_Debye_11;
				RE_index_X(i,j,k) = RE_Debye_11;
				RF_index_X(i,j,k) = RF_Debye_11;
				sigma_s_index_X(i,j,k) = sigma_Debye_11_S;
                
             elseif ( Brain_p(i,j,k) == 83)
				
                CB_index_X(i,j,k) = CB_Debye_12;

				RA_index_X_y(i,j,k) = RA_Debye_12;
				RA_index_X_z(i,j,k) = RA_Debye_12;

				RB_index_X_y(i,j,k) = RB_Debye_12;
				RB_index_X_z(i,j,k) = RB_Debye_12;

				RE_index_X_y(i,j,k) = RE_Debye_12;
				RE_index_X_z(i,j,k) = RE_Debye_12;

				RF_index_X_y(i,j,k) = RF_Debye_12;
				RF_index_X_z(i,j,k) = RF_Debye_12;

				RA_index_X(i,j,k) = RA_Debye_12;
				RB_index_X(i,j,k) = RB_Debye_12;
				RE_index_X(i,j,k) = RE_Debye_12;
				RF_index_X(i,j,k) = RF_Debye_12;
				sigma_s_index_X(i,j,k) = sigma_Debye_12_S;   
                
             elseif ( Brain_p(i,j,k) == 89 || Brain_p(i,j,k) == 95 || Brain_p(i,j,k) == 96 || Brain_p(i,j,k) == 101 || Brain_p(i,j,k) == 103 || Brain_p(i,j,k) == 105 || Brain_p(i,j,k) == 107 ...
                   || Brain_p(i,j,k) == 108 || Brain_p(i,j,k) == 109 || Brain_p(i,j,k) == 112 || Brain_p(i,j,k) == 114 || Brain_p(i,j,k) == 117 || Brain_p(i,j,k) == 118 || Brain_p(i,j,k) == 120 || Brain_p(i,j,k) == 124)
				
                CB_index_X(i,j,k) = CB_Debye_13;

				RA_index_X_y(i,j,k) = RA_Debye_13;
				RA_index_X_z(i,j,k) = RA_Debye_13;

				RB_index_X_y(i,j,k) = RB_Debye_13;
				RB_index_X_z(i,j,k) = RB_Debye_13;

				RE_index_X_y(i,j,k) = RE_Debye_13;
				RE_index_X_z(i,j,k) = RE_Debye_13;

				RF_index_X_y(i,j,k) = RF_Debye_13;
				RF_index_X_z(i,j,k) = RF_Debye_13;

				RA_index_X(i,j,k) = RA_Debye_13;
				RB_index_X(i,j,k) = RB_Debye_13;
				RE_index_X(i,j,k) = RE_Debye_13;
				RF_index_X(i,j,k) = RF_Debye_13;
				sigma_s_index_X(i,j,k) = sigma_Debye_13_S;  
                
             elseif ( Brain_p(i,j,k) == 26)
				
                CB_index_X(i,j,k) = CB_Debye_14;

				RA_index_X_y(i,j,k) = RA_Debye_14;
				RA_index_X_z(i,j,k) = RA_Debye_14;

				RB_index_X_y(i,j,k) = RB_Debye_14;
				RB_index_X_z(i,j,k) = RB_Debye_14;

				RE_index_X_y(i,j,k) = RE_Debye_14;
				RE_index_X_z(i,j,k) = RE_Debye_14;

				RF_index_X_y(i,j,k) = RF_Debye_14;
				RF_index_X_z(i,j,k) = RF_Debye_14;

				RA_index_X(i,j,k) = RA_Debye_14;
				RB_index_X(i,j,k) = RB_Debye_14;
				RE_index_X(i,j,k) = RE_Debye_14;
				RF_index_X(i,j,k) = RF_Debye_14;
				sigma_s_index_X(i,j,k) = sigma_Debye_14_S;   
                
             elseif ( Brain_p(i,j,k) == 97)
                 
				CB_index_X(i,j,k) = CB_Debye_15;

				RA_index_X_y(i,j,k) = RA_Debye_15;
				RA_index_X_z(i,j,k) = RA_Debye_15;

				RB_index_X_y(i,j,k) = RB_Debye_15;
				RB_index_X_z(i,j,k) = RB_Debye_15;

				RE_index_X_y(i,j,k) = RE_Debye_15;
				RE_index_X_z(i,j,k) = RE_Debye_15;

				RF_index_X_y(i,j,k) = RF_Debye_15;
				RF_index_X_z(i,j,k) = RF_Debye_15;

				RA_index_X(i,j,k) = RA_Debye_15;
				RB_index_X(i,j,k) = RB_Debye_15;
				RE_index_X(i,j,k) = RE_Debye_15;
				RF_index_X(i,j,k) = RF_Debye_15;
				sigma_s_index_X(i,j,k) = sigma_Debye_15_S;     
             
             elseif ( Brain_p(i,j,k) == 22 || Brain_p(i,j,k) == 98 || Brain_p(i,j,k) == 116)
				
                CB_index_X(i,j,k) = CB_Debye_16;

				RA_index_X_y(i,j,k) = RA_Debye_16;
				RA_index_X_z(i,j,k) = RA_Debye_16;

				RB_index_X_y(i,j,k) = RB_Debye_16;
				RB_index_X_z(i,j,k) = RB_Debye_16;

				RE_index_X_y(i,j,k) = RE_Debye_16;
				RE_index_X_z(i,j,k) = RE_Debye_16;

				RF_index_X_y(i,j,k) = RF_Debye_16;
				RF_index_X_z(i,j,k) = RF_Debye_16;

				RA_index_X(i,j,k) = RA_Debye_16;
				RB_index_X(i,j,k) = RB_Debye_16;
				RE_index_X(i,j,k) = RE_Debye_16;
				RF_index_X(i,j,k) = RF_Debye_16;
				sigma_s_index_X(i,j,k) = sigma_Debye_16_S;    
                
             elseif ( Brain_p(i,j,k) == 23)
				
                CB_index_X(i,j,k) = CB_Debye_17;

				RA_index_X_y(i,j,k) = RA_Debye_17;
				RA_index_X_z(i,j,k) = RA_Debye_17;

				RB_index_X_y(i,j,k) = RB_Debye_17;
				RB_index_X_z(i,j,k) = RB_Debye_17;

				RE_index_X_y(i,j,k) = RE_Debye_17;
				RE_index_X_z(i,j,k) = RE_Debye_17;

				RF_index_X_y(i,j,k) = RF_Debye_17;
				RF_index_X_z(i,j,k) = RF_Debye_17;

				RA_index_X(i,j,k) = RA_Debye_17;
				RB_index_X(i,j,k) = RB_Debye_17;
				RE_index_X(i,j,k) = RE_Debye_17;
				RF_index_X(i,j,k) = RF_Debye_17;
				sigma_s_index_X(i,j,k) = sigma_Debye_17_S;  
                
             elseif ( Brain_p(i,j,k) == 113)
				
                CB_index_X(i,j,k) = CB_Debye_18;

				RA_index_X_y(i,j,k) = RA_Debye_18;
				RA_index_X_z(i,j,k) = RA_Debye_18;

				RB_index_X_y(i,j,k) = RB_Debye_18;
				RB_index_X_z(i,j,k) = RB_Debye_18;

				RE_index_X_y(i,j,k) = RE_Debye_18;
				RE_index_X_z(i,j,k) = RE_Debye_18;

				RF_index_X_y(i,j,k) = RF_Debye_18;
				RF_index_X_z(i,j,k) = RF_Debye_18;

				RA_index_X(i,j,k) = RA_Debye_18;
				RB_index_X(i,j,k) = RB_Debye_18;
				RE_index_X(i,j,k) = RE_Debye_18;
				RF_index_X(i,j,k) = RF_Debye_18;
				sigma_s_index_X(i,j,k) = sigma_Debye_18_S;   
                
             elseif ( Brain_p(i,j,k) == 119)
				
                CB_index_X(i,j,k) = CB_Debye_19;

				RA_index_X_y(i,j,k) = RA_Debye_19;
				RA_index_X_z(i,j,k) = RA_Debye_19;

				RB_index_X_y(i,j,k) = RB_Debye_19;
				RB_index_X_z(i,j,k) = RB_Debye_19;

				RE_index_X_y(i,j,k) = RE_Debye_19;
				RE_index_X_z(i,j,k) = RE_Debye_19;

				RF_index_X_y(i,j,k) = RF_Debye_19;
				RF_index_X_z(i,j,k) = RF_Debye_19;

				RA_index_X(i,j,k) = RA_Debye_19;
				RB_index_X(i,j,k) = RB_Debye_19;
				RE_index_X(i,j,k) = RE_Debye_19;
				RF_index_X(i,j,k) = RF_Debye_19;
				sigma_s_index_X(i,j,k) = sigma_Debye_19_S;  
                
             elseif ( Brain_p(i,j,k) == 121)
				
                CB_index_X(i,j,k) = CB_Debye_20;

				RA_index_X_y(i,j,k) = RA_Debye_20;
				RA_index_X_z(i,j,k) = RA_Debye_20;

				RB_index_X_y(i,j,k) = RB_Debye_20;
				RB_index_X_z(i,j,k) = RB_Debye_20;

				RE_index_X_y(i,j,k) = RE_Debye_20;
				RE_index_X_z(i,j,k) = RE_Debye_20;

				RF_index_X_y(i,j,k) = RF_Debye_20;
				RF_index_X_z(i,j,k) = RF_Debye_20;

				RA_index_X(i,j,k) = RA_Debye_20;
				RB_index_X(i,j,k) = RB_Debye_20;
				RE_index_X(i,j,k) = RE_Debye_20;
				RF_index_X(i,j,k) = RF_Debye_20;
				sigma_s_index_X(i,j,k) = sigma_Debye_20_S;         
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
            
            if ( Brain_p(i,j,k) == 1)
                
				CB_index_Y(i,j,k) = CB_Debye_1;

				RA_index_Y_x(i,j,k) = RA_Debye_1;
				RA_index_Y_z(i,j,k) = RA_Debye_1;

				RB_index_Y_x(i,j,k) = RB_Debye_1;
				RB_index_Y_z(i,j,k) = RB_Debye_1;

				RE_index_Y_x(i,j,k) = RE_Debye_1;
				RE_index_Y_z(i,j,k) = RE_Debye_1;

				RF_index_Y_x(i,j,k) = RF_Debye_1;
				RF_index_Y_z(i,j,k) = RF_Debye_1;

				RA_index_Y(i,j,k) = RA_Debye_1;
				RB_index_Y(i,j,k) = RB_Debye_1;
				RE_index_Y(i,j,k) = RE_Debye_1;
				RF_index_Y(i,j,k) = RF_Debye_1;
				sigma_s_index_Y(i,j,k) = sigma_Debye_1_S;
                
            elseif ( Brain_p(i,j,k) == 2 || Brain_p(i,j,k) == 75 || Brain_p(i,j,k) == 92 || Brain_p(i,j,k) == 115 || Brain_p(i,j,k) == 122 || Brain_p(i,j,k) == 123)
               
                CB_index_Y(i,j,k) = CB_Debye_2;

				RA_index_Y_x(i,j,k) = RA_Debye_2;
				RA_index_Y_z(i,j,k) = RA_Debye_2;

				RB_index_Y_x(i,j,k) = RB_Debye_2;
				RB_index_Y_z(i,j,k) = RB_Debye_2;

				RE_index_Y_x(i,j,k) = RE_Debye_2;
				RE_index_Y_z(i,j,k) = RE_Debye_2;

				RF_index_Y_x(i,j,k) = RF_Debye_2;
				RF_index_Y_z(i,j,k) = RF_Debye_2;

				RA_index_Y(i,j,k) = RA_Debye_2;
				RB_index_Y(i,j,k) = RB_Debye_2;
				RE_index_Y(i,j,k) = RE_Debye_2;
				RF_index_Y(i,j,k) = RF_Debye_2;
				sigma_s_index_Y(i,j,k) = sigma_Debye_2_S;
                
             elseif ( Brain_p(i,j,k) == 3 || Brain_p(i,j,k) == 106 || Brain_p(i,j,k) == 111)
				
                CB_index_Y(i,j,k) = CB_Debye_3;

				RA_index_Y_x(i,j,k) = RA_Debye_3;
				RA_index_Y_z(i,j,k) = RA_Debye_3;

				RB_index_Y_x(i,j,k) = RB_Debye_3;
				RB_index_Y_z(i,j,k) = RB_Debye_3;

				RE_index_Y_x(i,j,k) = RE_Debye_3;
				RE_index_Y_z(i,j,k) = RE_Debye_3;

				RF_index_Y_x(i,j,k) = RF_Debye_3;
				RF_index_Y_z(i,j,k) = RF_Debye_3;

				RA_index_Y(i,j,k) = RA_Debye_3;
				RB_index_Y(i,j,k) = RB_Debye_3;
				RE_index_Y(i,j,k) = RE_Debye_3;
				RF_index_Y(i,j,k) = RF_Debye_3;
				sigma_s_index_Y(i,j,k) = sigma_Debye_3_S;
                
             elseif ( Brain_p(i,j,k) == 4 || Brain_p(i,j,k) == 5 || Brain_p(i,j,k) == 70 || Brain_p(i,j,k) == 71 || Brain_p(i,j,k) == 76 || Brain_p(i,j,k) == 81 || Brain_p(i,j,k) == 99 || Brain_p(i,j,k) == 125)
			
                CB_index_Y(i,j,k) = CB_Debye_4;

				RA_index_Y_x(i,j,k) = RA_Debye_4;
				RA_index_Y_z(i,j,k) = RA_Debye_4;

				RB_index_Y_x(i,j,k) = RB_Debye_4;
				RB_index_Y_z(i,j,k) = RB_Debye_4;

				RE_index_Y_x(i,j,k) = RE_Debye_4;
				RE_index_Y_z(i,j,k) = RE_Debye_4;

				RF_index_Y_x(i,j,k) = RF_Debye_4;
				RF_index_Y_z(i,j,k) = RF_Debye_4;

				RA_index_Y(i,j,k) = RA_Debye_4;
				RB_index_Y(i,j,k) = RB_Debye_4;
				RE_index_Y(i,j,k) = RE_Debye_4;
				RF_index_Y(i,j,k) = RF_Debye_4;
				sigma_s_index_Y(i,j,k) = sigma_Debye_4_S;
                
             elseif ( Brain_p(i,j,k) == 72)
				
                CB_index_Y(i,j,k) = CB_Debye_5;

				RA_index_Y_x(i,j,k) = RA_Debye_5;
				RA_index_Y_z(i,j,k) = RA_Debye_5;

				RB_index_Y_x(i,j,k) = RB_Debye_5;
				RB_index_Y_z(i,j,k) = RB_Debye_5;

				RE_index_Y_x(i,j,k) = RE_Debye_5;
				RE_index_Y_z(i,j,k) = RE_Debye_5;

				RF_index_Y_x(i,j,k) = RF_Debye_5;
				RF_index_Y_z(i,j,k) = RF_Debye_5;

				RA_index_Y(i,j,k) = RA_Debye_5;
				RB_index_Y(i,j,k) = RB_Debye_5;
				RE_index_Y(i,j,k) = RE_Debye_5;
				RF_index_Y(i,j,k) = RF_Debye_5;
				sigma_s_index_Y(i,j,k) = sigma_Debye_5_S;
                  
             elseif ( Brain_p(i,j,k) == 9 || Brain_p(i,j,k) == 102)
				
                CB_index_Y(i,j,k) = CB_Debye_6;

				RA_index_Y_x(i,j,k) = RA_Debye_6;
				RA_index_Y_z(i,j,k) = RA_Debye_6;

				RB_index_Y_x(i,j,k) = RB_Debye_6;
				RB_index_Y_z(i,j,k) = RB_Debye_6;

				RE_index_Y_x(i,j,k) = RE_Debye_6;
				RE_index_Y_z(i,j,k) = RE_Debye_6;

				RF_index_Y_x(i,j,k) = RF_Debye_6;
				RF_index_Y_z(i,j,k) = RF_Debye_6;

				RA_index_Y(i,j,k) = RA_Debye_6;
				RB_index_Y(i,j,k) = RB_Debye_6;
				RE_index_Y(i,j,k) = RE_Debye_6;
				RF_index_Y(i,j,k) = RF_Debye_6;
				sigma_s_index_Y(i,j,k) = sigma_Debye_6_S;
                
             elseif ( Brain_p(i,j,k) == 74)
				
                CB_index_Y(i,j,k) = CB_Debye_7;

				RA_index_Y_x(i,j,k) = RA_Debye_7;
				RA_index_Y_z(i,j,k) = RA_Debye_7;

				RB_index_Y_x(i,j,k) = RB_Debye_7;
				RB_index_Y_z(i,j,k) = RB_Debye_7;

				RE_index_Y_x(i,j,k) = RE_Debye_7;
				RE_index_Y_z(i,j,k) = RE_Debye_7;

				RF_index_Y_x(i,j,k) = RF_Debye_7;
				RF_index_Y_z(i,j,k) = RF_Debye_7;

				RA_index_Y(i,j,k) = RA_Debye_7;
				RB_index_Y(i,j,k) = RB_Debye_7;
				RE_index_Y(i,j,k) = RE_Debye_7;
				RF_index_Y(i,j,k) = RF_Debye_7;
				sigma_s_index_Y(i,j,k) = sigma_Debye_7_S;
                
             elseif ( Brain_p(i,j,k) == 77 || Brain_p(i,j,k) == 85 || Brain_p(i,j,k) == 91)
				
                CB_index_Y(i,j,k) = CB_Debye_8;

				RA_index_Y_x(i,j,k) = RA_Debye_8;
				RA_index_Y_z(i,j,k) = RA_Debye_8;

				RB_index_Y_x(i,j,k) = RB_Debye_8;
				RB_index_Y_z(i,j,k) = RB_Debye_8;

				RE_index_Y_x(i,j,k) = RE_Debye_8;
				RE_index_Y_z(i,j,k) = RE_Debye_8;

				RF_index_Y_x(i,j,k) = RF_Debye_8;
				RF_index_Y_z(i,j,k) = RF_Debye_8;

				RA_index_Y(i,j,k) = RA_Debye_8;
				RB_index_Y(i,j,k) = RB_Debye_8;
				RE_index_Y(i,j,k) = RE_Debye_8;
				RF_index_Y(i,j,k) = RF_Debye_8;
				sigma_s_index_Y(i,j,k) = sigma_Debye_8_S;  
                
             elseif ( Brain_p(i,j,k) == 78 || Brain_p(i,j,k) == 110)
				
                CB_index_Y(i,j,k) = CB_Debye_9;

				RA_index_Y_x(i,j,k) = RA_Debye_9;
				RA_index_Y_z(i,j,k) = RA_Debye_9;

				RB_index_Y_x(i,j,k) = RB_Debye_9;
				RB_index_Y_z(i,j,k) = RB_Debye_9;

				RE_index_Y_x(i,j,k) = RE_Debye_9;
				RE_index_Y_z(i,j,k) = RE_Debye_9;

				RF_index_Y_x(i,j,k) = RF_Debye_9;
				RF_index_Y_z(i,j,k) = RF_Debye_9;

				RA_index_Y(i,j,k) = RA_Debye_9;
				RB_index_Y(i,j,k) = RB_Debye_9;
				RE_index_Y(i,j,k) = RE_Debye_9;
				RF_index_Y(i,j,k) = RF_Debye_9;
				sigma_s_index_Y(i,j,k) = sigma_Debye_9_S;
                               
             elseif ( Brain_p(i,j,k) == 82 || Brain_p(i,j,k) == 30 || Brain_p(i,j,k) == 100)
				
                CB_index_Y(i,j,k) = CB_Debye_11;

				RA_index_Y_x(i,j,k) = RA_Debye_11;
				RA_index_Y_z(i,j,k) = RA_Debye_11;

				RB_index_Y_x(i,j,k) = RB_Debye_11;
				RB_index_Y_z(i,j,k) = RB_Debye_11;

				RE_index_Y_x(i,j,k) = RE_Debye_11;
				RE_index_Y_z(i,j,k) = RE_Debye_11;

				RF_index_Y_x(i,j,k) = RF_Debye_11;
				RF_index_Y_z(i,j,k) = RF_Debye_11;

				RA_index_Y(i,j,k) = RA_Debye_11;
				RB_index_Y(i,j,k) = RB_Debye_11;
				RE_index_Y(i,j,k) = RE_Debye_11;
				RF_index_Y(i,j,k) = RF_Debye_11;
				sigma_s_index_Y(i,j,k) = sigma_Debye_11_S;
                
             elseif ( Brain_p(i,j,k) == 83)
				
                CB_index_Y(i,j,k) = CB_Debye_12;

				RA_index_Y_x(i,j,k) = RA_Debye_12;
				RA_index_Y_z(i,j,k) = RA_Debye_12;

				RB_index_Y_x(i,j,k) = RB_Debye_12;
				RB_index_Y_z(i,j,k) = RB_Debye_12;

				RE_index_Y_x(i,j,k) = RE_Debye_12;
				RE_index_Y_z(i,j,k) = RE_Debye_12;

				RF_index_Y_x(i,j,k) = RF_Debye_12;
				RF_index_Y_z(i,j,k) = RF_Debye_12;

				RA_index_Y(i,j,k) = RA_Debye_12;
				RB_index_Y(i,j,k) = RB_Debye_12;
				RE_index_Y(i,j,k) = RE_Debye_12;
				RF_index_Y(i,j,k) = RF_Debye_12;
				sigma_s_index_Y(i,j,k) = sigma_Debye_12_S;   
                
             elseif ( Brain_p(i,j,k) == 89 || Brain_p(i,j,k) == 95 || Brain_p(i,j,k) == 96 || Brain_p(i,j,k) == 101 || Brain_p(i,j,k) == 103 || Brain_p(i,j,k) == 105 || Brain_p(i,j,k) == 107 ...
                   || Brain_p(i,j,k) == 108 || Brain_p(i,j,k) == 109 || Brain_p(i,j,k) == 112 || Brain_p(i,j,k) == 114 || Brain_p(i,j,k) == 117 || Brain_p(i,j,k) == 118 || Brain_p(i,j,k) == 120 || Brain_p(i,j,k) == 124)
				
                CB_index_Y(i,j,k) = CB_Debye_13;

				RA_index_Y_x(i,j,k) = RA_Debye_13;
				RA_index_Y_z(i,j,k) = RA_Debye_13;

				RB_index_Y_x(i,j,k) = RB_Debye_13;
				RB_index_Y_z(i,j,k) = RB_Debye_13;

				RE_index_Y_x(i,j,k) = RE_Debye_13;
				RE_index_Y_z(i,j,k) = RE_Debye_13;

				RF_index_Y_x(i,j,k) = RF_Debye_13;
				RF_index_Y_z(i,j,k) = RF_Debye_13;

				RA_index_Y(i,j,k) = RA_Debye_13;
				RB_index_Y(i,j,k) = RB_Debye_13;
				RE_index_Y(i,j,k) = RE_Debye_13;
				RF_index_Y(i,j,k) = RF_Debye_13;
				sigma_s_index_Y(i,j,k) = sigma_Debye_13_S;  
                
             elseif ( Brain_p(i,j,k) == 26)
				
                CB_index_Y(i,j,k) = CB_Debye_14;

				RA_index_Y_x(i,j,k) = RA_Debye_14;
				RA_index_Y_z(i,j,k) = RA_Debye_14;

				RB_index_Y_x(i,j,k) = RB_Debye_14;
				RB_index_Y_z(i,j,k) = RB_Debye_14;

				RE_index_Y_x(i,j,k) = RE_Debye_14;
				RE_index_Y_z(i,j,k) = RE_Debye_14;

				RF_index_Y_x(i,j,k) = RF_Debye_14;
				RF_index_Y_z(i,j,k) = RF_Debye_14;

				RA_index_Y(i,j,k) = RA_Debye_14;
				RB_index_Y(i,j,k) = RB_Debye_14;
				RE_index_Y(i,j,k) = RE_Debye_14;
				RF_index_Y(i,j,k) = RF_Debye_14;
				sigma_s_index_Y(i,j,k) = sigma_Debye_14_S;   
                
             elseif ( Brain_p(i,j,k) == 97)
                 
				CB_index_Y(i,j,k) = CB_Debye_15;

				RA_index_Y_x(i,j,k) = RA_Debye_15;
				RA_index_Y_z(i,j,k) = RA_Debye_15;

				RB_index_Y_x(i,j,k) = RB_Debye_15;
				RB_index_Y_z(i,j,k) = RB_Debye_15;

				RE_index_Y_x(i,j,k) = RE_Debye_15;
				RE_index_Y_z(i,j,k) = RE_Debye_15;

				RF_index_Y_x(i,j,k) = RF_Debye_15;
				RF_index_Y_z(i,j,k) = RF_Debye_15;

				RA_index_Y(i,j,k) = RA_Debye_15;
				RB_index_Y(i,j,k) = RB_Debye_15;
				RE_index_Y(i,j,k) = RE_Debye_15;
				RF_index_Y(i,j,k) = RF_Debye_15;
				sigma_s_index_Y(i,j,k) = sigma_Debye_15_S;     
             
             elseif ( Brain_p(i,j,k) == 22 || Brain_p(i,j,k) == 98 || Brain_p(i,j,k) == 116)
				
                CB_index_Y(i,j,k) = CB_Debye_16;

				RA_index_Y_x(i,j,k) = RA_Debye_16;
				RA_index_Y_z(i,j,k) = RA_Debye_16;

				RB_index_Y_x(i,j,k) = RB_Debye_16;
				RB_index_Y_z(i,j,k) = RB_Debye_16;

				RE_index_Y_x(i,j,k) = RE_Debye_16;
				RE_index_Y_z(i,j,k) = RE_Debye_16;

				RF_index_Y_x(i,j,k) = RF_Debye_16;
				RF_index_Y_z(i,j,k) = RF_Debye_16;

				RA_index_Y(i,j,k) = RA_Debye_16;
				RB_index_Y(i,j,k) = RB_Debye_16;
				RE_index_Y(i,j,k) = RE_Debye_16;
				RF_index_Y(i,j,k) = RF_Debye_16;
				sigma_s_index_Y(i,j,k) = sigma_Debye_16_S;    
                
             elseif ( Brain_p(i,j,k) == 23)
				
                CB_index_Y(i,j,k) = CB_Debye_17;

				RA_index_Y_x(i,j,k) = RA_Debye_17;
				RA_index_Y_z(i,j,k) = RA_Debye_17;

				RB_index_Y_x(i,j,k) = RB_Debye_17;
				RB_index_Y_z(i,j,k) = RB_Debye_17;

				RE_index_Y_x(i,j,k) = RE_Debye_17;
				RE_index_Y_z(i,j,k) = RE_Debye_17;

				RF_index_Y_x(i,j,k) = RF_Debye_17;
				RF_index_Y_z(i,j,k) = RF_Debye_17;

				RA_index_Y(i,j,k) = RA_Debye_17;
				RB_index_Y(i,j,k) = RB_Debye_17;
				RE_index_Y(i,j,k) = RE_Debye_17;
				RF_index_Y(i,j,k) = RF_Debye_17;
				sigma_s_index_Y(i,j,k) = sigma_Debye_17_S;  
                
             elseif ( Brain_p(i,j,k) == 113)
				
                CB_index_Y(i,j,k) = CB_Debye_18;

				RA_index_Y_x(i,j,k) = RA_Debye_18;
				RA_index_Y_z(i,j,k) = RA_Debye_18;

				RB_index_Y_x(i,j,k) = RB_Debye_18;
				RB_index_Y_z(i,j,k) = RB_Debye_18;

				RE_index_Y_x(i,j,k) = RE_Debye_18;
				RE_index_Y_z(i,j,k) = RE_Debye_18;

				RF_index_Y_x(i,j,k) = RF_Debye_18;
				RF_index_Y_z(i,j,k) = RF_Debye_18;

				RA_index_Y(i,j,k) = RA_Debye_18;
				RB_index_Y(i,j,k) = RB_Debye_18;
				RE_index_Y(i,j,k) = RE_Debye_18;
				RF_index_Y(i,j,k) = RF_Debye_18;
				sigma_s_index_Y(i,j,k) = sigma_Debye_18_S;   
                
             elseif ( Brain_p(i,j,k) == 119)
				
                CB_index_Y(i,j,k) = CB_Debye_19;

				RA_index_Y_x(i,j,k) = RA_Debye_19;
				RA_index_Y_z(i,j,k) = RA_Debye_19;

				RB_index_Y_x(i,j,k) = RB_Debye_19;
				RB_index_Y_z(i,j,k) = RB_Debye_19;

				RE_index_Y_x(i,j,k) = RE_Debye_19;
				RE_index_Y_z(i,j,k) = RE_Debye_19;

				RF_index_Y_x(i,j,k) = RF_Debye_19;
				RF_index_Y_z(i,j,k) = RF_Debye_19;

				RA_index_Y(i,j,k) = RA_Debye_19;
				RB_index_Y(i,j,k) = RB_Debye_19;
				RE_index_Y(i,j,k) = RE_Debye_19;
				RF_index_Y(i,j,k) = RF_Debye_19;
				sigma_s_index_Y(i,j,k) = sigma_Debye_19_S;  
                
             elseif ( Brain_p(i,j,k) == 121)
				
                CB_index_Y(i,j,k) = CB_Debye_20;

				RA_index_Y_x(i,j,k) = RA_Debye_20;
				RA_index_Y_z(i,j,k) = RA_Debye_20;

				RB_index_Y_x(i,j,k) = RB_Debye_20;
				RB_index_Y_z(i,j,k) = RB_Debye_20;

				RE_index_Y_x(i,j,k) = RE_Debye_20;
				RE_index_Y_z(i,j,k) = RE_Debye_20;

				RF_index_Y_x(i,j,k) = RF_Debye_20;
				RF_index_Y_z(i,j,k) = RF_Debye_20;

				RA_index_Y(i,j,k) = RA_Debye_20;
				RB_index_Y(i,j,k) = RB_Debye_20;
				RE_index_Y(i,j,k) = RE_Debye_20;
				RF_index_Y(i,j,k) = RF_Debye_20;
				sigma_s_index_Y(i,j,k) = sigma_Debye_20_S;         
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
       
              
           if ( Brain_p(i,j,k) == 1)
                
				CB_index_Z(i,j,k) = CB_Debye_1;

				RA_index_Z_x(i,j,k) = RA_Debye_1;
				RA_index_Z_y(i,j,k) = RA_Debye_1;

				RB_index_Z_x(i,j,k) = RB_Debye_1;
				RB_index_Z_y(i,j,k) = RB_Debye_1;

				RE_index_Z_x(i,j,k) = RE_Debye_1;
				RE_index_Z_y(i,j,k) = RE_Debye_1;

				RF_index_Z_x(i,j,k) = RF_Debye_1;
				RF_index_Z_y(i,j,k) = RF_Debye_1;

				RA_index_Z(i,j,k) = RA_Debye_1;
				RB_index_Z(i,j,k) = RB_Debye_1;
				RE_index_Z(i,j,k) = RE_Debye_1;
				RF_index_Z(i,j,k) = RF_Debye_1;
				sigma_s_index_Z(i,j,k) = sigma_Debye_1_S;
                
            elseif ( Brain_p(i,j,k) == 2 || Brain_p(i,j,k) == 75 || Brain_p(i,j,k) == 92 || Brain_p(i,j,k) == 115 || Brain_p(i,j,k) == 122 || Brain_p(i,j,k) == 123)
               
                CB_index_Z(i,j,k) = CB_Debye_2;

				RA_index_Z_x(i,j,k) = RA_Debye_2;
				RA_index_Z_y(i,j,k) = RA_Debye_2;

				RB_index_Z_x(i,j,k) = RB_Debye_2;
				RB_index_Z_y(i,j,k) = RB_Debye_2;

				RE_index_Z_x(i,j,k) = RE_Debye_2;
				RE_index_Z_y(i,j,k) = RE_Debye_2;

				RF_index_Z_x(i,j,k) = RF_Debye_2;
				RF_index_Z_y(i,j,k) = RF_Debye_2;

				RA_index_Z(i,j,k) = RA_Debye_2;
				RB_index_Z(i,j,k) = RB_Debye_2;
				RE_index_Z(i,j,k) = RE_Debye_2;
				RF_index_Z(i,j,k) = RF_Debye_2;
				sigma_s_index_Z(i,j,k) = sigma_Debye_2_S;
                
             elseif ( Brain_p(i,j,k) == 3 || Brain_p(i,j,k) == 106 || Brain_p(i,j,k) == 111)
				
                CB_index_Z(i,j,k) = CB_Debye_3;

				RA_index_Z_x(i,j,k) = RA_Debye_3;
				RA_index_Z_y(i,j,k) = RA_Debye_3;

				RB_index_Z_x(i,j,k) = RB_Debye_3;
				RB_index_Z_y(i,j,k) = RB_Debye_3;

				RE_index_Z_x(i,j,k) = RE_Debye_3;
				RE_index_Z_y(i,j,k) = RE_Debye_3;

				RF_index_Z_x(i,j,k) = RF_Debye_3;
				RF_index_Z_y(i,j,k) = RF_Debye_3;

				RA_index_Z(i,j,k) = RA_Debye_3;
				RB_index_Z(i,j,k) = RB_Debye_3;
				RE_index_Z(i,j,k) = RE_Debye_3;
				RF_index_Z(i,j,k) = RF_Debye_3;
				sigma_s_index_Z(i,j,k) = sigma_Debye_3_S;
                
             elseif ( Brain_p(i,j,k) == 4 || Brain_p(i,j,k) == 5 || Brain_p(i,j,k) == 70 || Brain_p(i,j,k) == 71 || Brain_p(i,j,k) == 76 || Brain_p(i,j,k) == 81 || Brain_p(i,j,k) == 99 || Brain_p(i,j,k) == 125)
			
                CB_index_Z(i,j,k) = CB_Debye_4;

				RA_index_Z_x(i,j,k) = RA_Debye_4;
				RA_index_Z_y(i,j,k) = RA_Debye_4;

				RB_index_Z_x(i,j,k) = RB_Debye_4;
				RB_index_Z_y(i,j,k) = RB_Debye_4;

				RE_index_Z_x(i,j,k) = RE_Debye_4;
				RE_index_Z_y(i,j,k) = RE_Debye_4;

				RF_index_Z_x(i,j,k) = RF_Debye_4;
				RF_index_Z_y(i,j,k) = RF_Debye_4;

				RA_index_Z(i,j,k) = RA_Debye_4;
				RB_index_Z(i,j,k) = RB_Debye_4;
				RE_index_Z(i,j,k) = RE_Debye_4;
				RF_index_Z(i,j,k) = RF_Debye_4;
				sigma_s_index_Z(i,j,k) = sigma_Debye_4_S;
                
             elseif ( Brain_p(i,j,k) == 72)
				
                CB_index_Z(i,j,k) = CB_Debye_5;

				RA_index_Z_x(i,j,k) = RA_Debye_5;
				RA_index_Z_y(i,j,k) = RA_Debye_5;

				RB_index_Z_x(i,j,k) = RB_Debye_5;
				RB_index_Z_y(i,j,k) = RB_Debye_5;

				RE_index_Z_x(i,j,k) = RE_Debye_5;
				RE_index_Z_y(i,j,k) = RE_Debye_5;

				RF_index_Z_x(i,j,k) = RF_Debye_5;
				RF_index_Z_y(i,j,k) = RF_Debye_5;

				RA_index_Z(i,j,k) = RA_Debye_5;
				RB_index_Z(i,j,k) = RB_Debye_5;
				RE_index_Z(i,j,k) = RE_Debye_5;
				RF_index_Z(i,j,k) = RF_Debye_5;
				sigma_s_index_Z(i,j,k) = sigma_Debye_5_S;
                  
             elseif ( Brain_p(i,j,k) == 9 || Brain_p(i,j,k) == 102)
				
                CB_index_Z(i,j,k) = CB_Debye_6;

				RA_index_Z_x(i,j,k) = RA_Debye_6;
				RA_index_Z_y(i,j,k) = RA_Debye_6;

				RB_index_Z_x(i,j,k) = RB_Debye_6;
				RB_index_Z_y(i,j,k) = RB_Debye_6;

				RE_index_Z_x(i,j,k) = RE_Debye_6;
				RE_index_Z_y(i,j,k) = RE_Debye_6;

				RF_index_Z_x(i,j,k) = RF_Debye_6;
				RF_index_Z_y(i,j,k) = RF_Debye_6;

				RA_index_Z(i,j,k) = RA_Debye_6;
				RB_index_Z(i,j,k) = RB_Debye_6;
				RE_index_Z(i,j,k) = RE_Debye_6;
				RF_index_Z(i,j,k) = RF_Debye_6;
				sigma_s_index_Z(i,j,k) = sigma_Debye_6_S;
                
             elseif ( Brain_p(i,j,k) == 74)
				
                CB_index_Z(i,j,k) = CB_Debye_7;

				RA_index_Z_x(i,j,k) = RA_Debye_7;
				RA_index_Z_y(i,j,k) = RA_Debye_7;

				RB_index_Z_x(i,j,k) = RB_Debye_7;
				RB_index_Z_y(i,j,k) = RB_Debye_7;

				RE_index_Z_x(i,j,k) = RE_Debye_7;
				RE_index_Z_y(i,j,k) = RE_Debye_7;

				RF_index_Z_x(i,j,k) = RF_Debye_7;
				RF_index_Z_y(i,j,k) = RF_Debye_7;

				RA_index_Z(i,j,k) = RA_Debye_7;
				RB_index_Z(i,j,k) = RB_Debye_7;
				RE_index_Z(i,j,k) = RE_Debye_7;
				RF_index_Z(i,j,k) = RF_Debye_7;
				sigma_s_index_Z(i,j,k) = sigma_Debye_7_S;
                
             elseif ( Brain_p(i,j,k) == 77 || Brain_p(i,j,k) == 85 || Brain_p(i,j,k) == 91)
				
                CB_index_Z(i,j,k) = CB_Debye_8;

				RA_index_Z_x(i,j,k) = RA_Debye_8;
				RA_index_Z_y(i,j,k) = RA_Debye_8;

				RB_index_Z_x(i,j,k) = RB_Debye_8;
				RB_index_Z_y(i,j,k) = RB_Debye_8;

				RE_index_Z_x(i,j,k) = RE_Debye_8;
				RE_index_Z_y(i,j,k) = RE_Debye_8;

				RF_index_Z_x(i,j,k) = RF_Debye_8;
				RF_index_Z_y(i,j,k) = RF_Debye_8;

				RA_index_Z(i,j,k) = RA_Debye_8;
				RB_index_Z(i,j,k) = RB_Debye_8;
				RE_index_Z(i,j,k) = RE_Debye_8;
				RF_index_Z(i,j,k) = RF_Debye_8;
				sigma_s_index_Z(i,j,k) = sigma_Debye_8_S;  
                
             elseif ( Brain_p(i,j,k) == 78 || Brain_p(i,j,k) == 110)
				
                CB_index_Z(i,j,k) = CB_Debye_9;

				RA_index_Z_x(i,j,k) = RA_Debye_9;
				RA_index_Z_y(i,j,k) = RA_Debye_9;

				RB_index_Z_x(i,j,k) = RB_Debye_9;
				RB_index_Z_y(i,j,k) = RB_Debye_9;

				RE_index_Z_x(i,j,k) = RE_Debye_9;
				RE_index_Z_y(i,j,k) = RE_Debye_9;

				RF_index_Z_x(i,j,k) = RF_Debye_9;
				RF_index_Z_y(i,j,k) = RF_Debye_9;

				RA_index_Z(i,j,k) = RA_Debye_9;
				RB_index_Z(i,j,k) = RB_Debye_9;
				RE_index_Z(i,j,k) = RE_Debye_9;
				RF_index_Z(i,j,k) = RF_Debye_9;
				sigma_s_index_Z(i,j,k) = sigma_Debye_9_S;
                                
             elseif ( Brain_p(i,j,k) == 82 || Brain_p(i,j,k) == 30 || Brain_p(i,j,k) == 100)
				
                CB_index_Z(i,j,k) = CB_Debye_11;

				RA_index_Z_x(i,j,k) = RA_Debye_11;
				RA_index_Z_y(i,j,k) = RA_Debye_11;

				RB_index_Z_x(i,j,k) = RB_Debye_11;
				RB_index_Z_y(i,j,k) = RB_Debye_11;

				RE_index_Z_x(i,j,k) = RE_Debye_11;
				RE_index_Z_y(i,j,k) = RE_Debye_11;

				RF_index_Z_x(i,j,k) = RF_Debye_11;
				RF_index_Z_y(i,j,k) = RF_Debye_11;

				RA_index_Z(i,j,k) = RA_Debye_11;
				RB_index_Z(i,j,k) = RB_Debye_11;
				RE_index_Z(i,j,k) = RE_Debye_11;
				RF_index_Z(i,j,k) = RF_Debye_11;
				sigma_s_index_Z(i,j,k) = sigma_Debye_11_S;
                
             elseif ( Brain_p(i,j,k) == 83)
				
                CB_index_Z(i,j,k) = CB_Debye_12;

				RA_index_Z_x(i,j,k) = RA_Debye_12;
				RA_index_Z_y(i,j,k) = RA_Debye_12;

				RB_index_Z_x(i,j,k) = RB_Debye_12;
				RB_index_Z_y(i,j,k) = RB_Debye_12;

				RE_index_Z_x(i,j,k) = RE_Debye_12;
				RE_index_Z_y(i,j,k) = RE_Debye_12;

				RF_index_Z_x(i,j,k) = RF_Debye_12;
				RF_index_Z_y(i,j,k) = RF_Debye_12;

				RA_index_Z(i,j,k) = RA_Debye_12;
				RB_index_Z(i,j,k) = RB_Debye_12;
				RE_index_Z(i,j,k) = RE_Debye_12;
				RF_index_Z(i,j,k) = RF_Debye_12;
				sigma_s_index_Z(i,j,k) = sigma_Debye_12_S;   
                
             elseif ( Brain_p(i,j,k) == 89 || Brain_p(i,j,k) == 95 || Brain_p(i,j,k) == 96 || Brain_p(i,j,k) == 101 || Brain_p(i,j,k) == 103 || Brain_p(i,j,k) == 105 || Brain_p(i,j,k) == 107 ...
                   || Brain_p(i,j,k) == 108 || Brain_p(i,j,k) == 109 || Brain_p(i,j,k) == 112 || Brain_p(i,j,k) == 114 || Brain_p(i,j,k) == 117 || Brain_p(i,j,k) == 118 || Brain_p(i,j,k) == 120 || Brain_p(i,j,k) == 124)
				
                CB_index_Z(i,j,k) = CB_Debye_13;

				RA_index_Z_x(i,j,k) = RA_Debye_13;
				RA_index_Z_y(i,j,k) = RA_Debye_13;

				RB_index_Z_x(i,j,k) = RB_Debye_13;
				RB_index_Z_y(i,j,k) = RB_Debye_13;

				RE_index_Z_x(i,j,k) = RE_Debye_13;
				RE_index_Z_y(i,j,k) = RE_Debye_13;

				RF_index_Z_x(i,j,k) = RF_Debye_13;
				RF_index_Z_y(i,j,k) = RF_Debye_13;

				RA_index_Z(i,j,k) = RA_Debye_13;
				RB_index_Z(i,j,k) = RB_Debye_13;
				RE_index_Z(i,j,k) = RE_Debye_13;
				RF_index_Z(i,j,k) = RF_Debye_13;
				sigma_s_index_Z(i,j,k) = sigma_Debye_13_S;  
                
             elseif ( Brain_p(i,j,k) == 26)
				
                CB_index_Z(i,j,k) = CB_Debye_14;

				RA_index_Z_x(i,j,k) = RA_Debye_14;
				RA_index_Z_y(i,j,k) = RA_Debye_14;

				RB_index_Z_x(i,j,k) = RB_Debye_14;
				RB_index_Z_y(i,j,k) = RB_Debye_14;

				RE_index_Z_x(i,j,k) = RE_Debye_14;
				RE_index_Z_y(i,j,k) = RE_Debye_14;

				RF_index_Z_x(i,j,k) = RF_Debye_14;
				RF_index_Z_y(i,j,k) = RF_Debye_14;

				RA_index_Z(i,j,k) = RA_Debye_14;
				RB_index_Z(i,j,k) = RB_Debye_14;
				RE_index_Z(i,j,k) = RE_Debye_14;
				RF_index_Z(i,j,k) = RF_Debye_14;
				sigma_s_index_Z(i,j,k) = sigma_Debye_14_S;   
                
             elseif ( Brain_p(i,j,k) == 97)
                 
				CB_index_Z(i,j,k) = CB_Debye_15;

				RA_index_Z_x(i,j,k) = RA_Debye_15;
				RA_index_Z_y(i,j,k) = RA_Debye_15;

				RB_index_Z_x(i,j,k) = RB_Debye_15;
				RB_index_Z_y(i,j,k) = RB_Debye_15;

				RE_index_Z_x(i,j,k) = RE_Debye_15;
				RE_index_Z_y(i,j,k) = RE_Debye_15;

				RF_index_Z_x(i,j,k) = RF_Debye_15;
				RF_index_Z_y(i,j,k) = RF_Debye_15;

				RA_index_Z(i,j,k) = RA_Debye_15;
				RB_index_Z(i,j,k) = RB_Debye_15;
				RE_index_Z(i,j,k) = RE_Debye_15;
				RF_index_Z(i,j,k) = RF_Debye_15;
				sigma_s_index_Z(i,j,k) = sigma_Debye_15_S;     
             
             elseif ( Brain_p(i,j,k) == 22 || Brain_p(i,j,k) == 98 || Brain_p(i,j,k) == 116)
				
                CB_index_Z(i,j,k) = CB_Debye_16;

				RA_index_Z_x(i,j,k) = RA_Debye_16;
				RA_index_Z_y(i,j,k) = RA_Debye_16;

				RB_index_Z_x(i,j,k) = RB_Debye_16;
				RB_index_Z_y(i,j,k) = RB_Debye_16;

				RE_index_Z_x(i,j,k) = RE_Debye_16;
				RE_index_Z_y(i,j,k) = RE_Debye_16;

				RF_index_Z_x(i,j,k) = RF_Debye_16;
				RF_index_Z_y(i,j,k) = RF_Debye_16;

				RA_index_Z(i,j,k) = RA_Debye_16;
				RB_index_Z(i,j,k) = RB_Debye_16;
				RE_index_Z(i,j,k) = RE_Debye_16;
				RF_index_Z(i,j,k) = RF_Debye_16;
				sigma_s_index_Z(i,j,k) = sigma_Debye_16_S;    
                
             elseif ( Brain_p(i,j,k) == 23)
				
                CB_index_Z(i,j,k) = CB_Debye_17;

				RA_index_Z_x(i,j,k) = RA_Debye_17;
				RA_index_Z_y(i,j,k) = RA_Debye_17;

				RB_index_Z_x(i,j,k) = RB_Debye_17;
				RB_index_Z_y(i,j,k) = RB_Debye_17;

				RE_index_Z_x(i,j,k) = RE_Debye_17;
				RE_index_Z_y(i,j,k) = RE_Debye_17;

				RF_index_Z_x(i,j,k) = RF_Debye_17;
				RF_index_Z_y(i,j,k) = RF_Debye_17;

				RA_index_Z(i,j,k) = RA_Debye_17;
				RB_index_Z(i,j,k) = RB_Debye_17;
				RE_index_Z(i,j,k) = RE_Debye_17;
				RF_index_Z(i,j,k) = RF_Debye_17;
				sigma_s_index_Z(i,j,k) = sigma_Debye_17_S;  
                
             elseif ( Brain_p(i,j,k) == 113)
				
                CB_index_Z(i,j,k) = CB_Debye_18;

				RA_index_Z_x(i,j,k) = RA_Debye_18;
				RA_index_Z_y(i,j,k) = RA_Debye_18;

				RB_index_Z_x(i,j,k) = RB_Debye_18;
				RB_index_Z_y(i,j,k) = RB_Debye_18;

				RE_index_Z_x(i,j,k) = RE_Debye_18;
				RE_index_Z_y(i,j,k) = RE_Debye_18;

				RF_index_Z_x(i,j,k) = RF_Debye_18;
				RF_index_Z_y(i,j,k) = RF_Debye_18;

				RA_index_Z(i,j,k) = RA_Debye_18;
				RB_index_Z(i,j,k) = RB_Debye_18;
				RE_index_Z(i,j,k) = RE_Debye_18;
				RF_index_Z(i,j,k) = RF_Debye_18;
				sigma_s_index_Z(i,j,k) = sigma_Debye_18_S;   
                
             elseif ( Brain_p(i,j,k) == 119)
				
                CB_index_Z(i,j,k) = CB_Debye_19;

				RA_index_Z_x(i,j,k) = RA_Debye_19;
				RA_index_Z_y(i,j,k) = RA_Debye_19;

				RB_index_Z_x(i,j,k) = RB_Debye_19;
				RB_index_Z_y(i,j,k) = RB_Debye_19;

				RE_index_Z_x(i,j,k) = RE_Debye_19;
				RE_index_Z_y(i,j,k) = RE_Debye_19;

				RF_index_Z_x(i,j,k) = RF_Debye_19;
				RF_index_Z_y(i,j,k) = RF_Debye_19;

				RA_index_Z(i,j,k) = RA_Debye_19;
				RB_index_Z(i,j,k) = RB_Debye_19;
				RE_index_Z(i,j,k) = RE_Debye_19;
				RF_index_Z(i,j,k) = RF_Debye_19;
				sigma_s_index_Z(i,j,k) = sigma_Debye_19_S;  
                
             elseif ( Brain_p(i,j,k) == 121)
				
                CB_index_Z(i,j,k) = CB_Debye_20;

				RA_index_Z_x(i,j,k) = RA_Debye_20;
				RA_index_Z_y(i,j,k) = RA_Debye_20;

				RB_index_Z_x(i,j,k) = RB_Debye_20;
				RB_index_Z_y(i,j,k) = RB_Debye_20;

				RE_index_Z_x(i,j,k) = RE_Debye_20;
				RE_index_Z_y(i,j,k) = RE_Debye_20;

				RF_index_Z_x(i,j,k) = RF_Debye_20;
				RF_index_Z_y(i,j,k) = RF_Debye_20;

				RA_index_Z(i,j,k) = RA_Debye_20;
				RB_index_Z(i,j,k) = RB_Debye_20;
				RE_index_Z(i,j,k) = RE_Debye_20;
				RF_index_Z(i,j,k) = RF_Debye_20;
				sigma_s_index_Z(i,j,k) = sigma_Debye_20_S;         
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

C = {'v1_Ex','v2_Ex','Hz','Hy','Ex_1','Ex',...
'CA_Tx','CB_Tx','CB_RA_X_y','CB_RA_X_z',...
'CB_RB_X_y','CB_RB_X_z','GEZ_Sxy','GEZ_Sxz','RE_index_X_y','RF_index_X_y',...
'RE_index_X_z','RF_index_X_z',...
'CB_RA_sigma_s_X_X', 'CB_RB_sigma_s_X_X', 'GP_Sx', 'RE_index_X', 'RF_index_X',...
'v1_Ey','v2_Ey','Hx','Ey_1','Ey',...
'CA_Ty','CB_Ty','CB_RA_Y_z','CB_RA_Y_x',...
'CB_RB_Y_z','CB_RB_Y_x','GEZ_Syz','GEZ_Syx','RE_index_Y_z','RF_index_Y_z',...
'RE_index_Y_x','RF_index_Y_x',...
'CB_RA_sigma_s_Y_Y', 'CB_RB_sigma_s_Y_Y', 'GP_Sy', 'RE_index_Y', 'RF_index_Y',...
'v1_Ez','v2_Ez','Ez_1','Ez',...
'CA_Tz','CB_Tz','CB_RA_Z_x','CB_RA_Z_y',...
'CB_RB_Z_x','CB_RB_Z_y','GEZ_Szx','GEZ_Szy','RE_index_Z_x','RF_index_Z_x',...
'RE_index_Z_y','RF_index_Z_y',...
'CB_RA_sigma_s_Z_Z', 'CB_RB_sigma_s_Z_Z', 'GP_Sz', 'RE_index_Z', 'RF_index_Z',...
'v1_Hx','v2_Hx','Ez','Ey','RAH_index_X_y','RBH_index_X_y','Hxy','RAH_index_X_z','RBH_index_X_z','Hxz',...
'REH_index_X_y','RFH_index_X_y','REH_index_X_z','RFH_index_X_z',...
'v1_Hy','v2_Hy','Ex','RAH_index_Y_z','RBH_index_Y_z','Hyz','RAH_index_Y_x','RBH_index_Y_x','Hyx',...
'REH_index_Y_z','RFH_index_Y_z','REH_index_Y_x','RFH_index_Y_x',...
'v1_Hz','v2_Hz','RAH_index_Z_x','RBH_index_Z_x','Hzx','RAH_index_Z_y','RBH_index_Z_y','Hzy',...
'REH_index_Z_x','RFH_index_Z_x','REH_index_Z_y','RFH_index_Z_y',};

for n=1:length(C)
    str = strcat(C{n},'=gpuArray(',C{n},');');
    eval(str);
end

%% ourside the cycle

% calculate the SAR and electric power for every time step and the whole
% simulation time (RI method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_power_single_RI_3 = zeros(Imax+1, Jmax+1, Kmax+1,'gpuArray');
Brain_SAR_single_RI_3 =zeros(Imax+1, Jmax+1, Kmax+1,'gpuArray');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_power_single_RI_4 = zeros(Imax+1, Jmax+1, Kmax+1,'gpuArray');
Brain_SAR_single_RI_4 =zeros(Imax+1, Jmax+1, Kmax+1,'gpuArray');
  
E_power_total_RI_4 = zeros(Imax+1, Jmax+1, Kmax+1,'gpuArray');
Brain_SAR_total_RI_4 = zeros(Imax+1, Jmax+1, Kmax+1,'gpuArray');


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
		  
    % Auxiliary variable
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

	% Auxiliary variable
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
    
	% Auxiliary variable
	GEZ_Szx(i,j,k) = RE_index_Z_x(i,j,k) .* GEZ_Szx(i,j,k) + RF_index_Z_x(i,j,k) .* v1_Ez(i,j,k);
	GEZ_Szy(i,j,k) = RE_index_Z_y(i,j,k) .* GEZ_Szy(i,j,k) + RF_index_Z_y(i,j,k) .* v2_Ez(i,j,k);
	GP_Sz(i,j,k) = RE_index_Z(i,j,k) .* GP_Sz(i,j,k) + RF_index_Z(i,j,k) .* Ez_1(i,j,k);

%%%%   源
%   source(n) = exp(-(((n + 1) * EM_DELTAT - t0) / (SPREAD))^2);
%   Ez(isource01+0,jsource01-0,ksource01)= Ez(isource01+0,jsource01-0,ksource01) + source(n);
source(n) = -2.0 * (((n + 1 - 1) * EM_DELTAT - T0) / (SPREAD)) * exp(-(((n + 1 - 1) * EM_DELTAT - T0) / (SPREAD))^2);

Ez(IC,JC) = Ez(IC,JC) + source(n);

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
  Hzy(i,j,k) = REH_index_Z_y(i, ...
      j,k) .* Hzy(i,j,k) + RFH_index_Z_y(i,j,k) .* v2_Hz(i,j,k);

%***********************************************************************
% get the SAR in the brain for each time step and the whole simlulation
% time (RI method)
%***********************************************************************
 %% Data post processing   
 
% Probe Point   
    XGD_EH_RIPML_Debye_RI_77(n) = Ez(IC,JC - 20,KC);
     XGD_EH_RIPML_Debye_RI_1_77(n) = Ez(IC,JC + 20,KC + 20);

%***********************************************************************
% get the SAR in the brain for each time step and the whole simlulation
% time (RI method)
%***********************************************************************

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if n<=6000
       
     i=1:Imax+1;
     j=1:Jmax+1;
     k=1:Kmax+1;
  
     E_power_single_RI_3(i,j,k)= Ex(i,j,k).^2 + Ey(i,j,k).^2 + Ez(i,j,k).^2;
     Brain_SAR_single_RI_3(i,j,k)=index_sig(i,j,k).* E_power_single_RI_3(i,j,k)./index_rou(i,j,k);
  
     
  end
  
  if n<=8000
       
     i=1:Imax+1;
     j=1:Jmax+1;
     k=1:Kmax+1;
  
     E_power_single_RI_4(i,j,k)= Ex(i,j,k).^2 + Ey(i,j,k).^2 + Ez(i,j,k).^2;
     Brain_SAR_single_RI_4(i,j,k)=index_sig(i,j,k).* E_power_single_RI_4(i,j,k)./index_rou(i,j,k);
  
     E_power_total_RI_4(i,j,k)= E_power_total_RI_4(i,j,k) + E_power_single_RI_4(i,j,k);
     Brain_SAR_total_RI_4(i,j,k)=Brain_SAR_total_RI_4(i,j,k) + Brain_SAR_single_RI_4(i,j,k);
     
  end
  

  %}

end
toc;

[XGD_EH_RIPML_Debye_RI_77,Brain_SAR_single_RI_3,Brain_SAR_single_RI_4] = gather(XGD_EH_RIPML_Debye_RI_77,Brain_SAR_single_RI_3,Brain_SAR_single_RI_4);
save('./RI.mat','XGD_EH_RIPML_Debye_RI_77','Brain_SAR_single_RI_3','Brain_SAR_single_RI_4');
[Brain_SAR_total_RI_4]=gather(Brain_SAR_total_RI_4);
save('./RI.mat','Brain_SAR_total_RI_4','-append');
[E_power_single_RI_3,E_power_single_RI_4]=gather(E_power_single_RI_3,E_power_single_RI_4);
save('./RI.mat','E_power_single_RI_3','E_power_single_RI_4','-append');
[E_power_total_RI_4]=gather(E_power_total_RI_4);
save('./RI.mat','E_power_total_RI_4','-append');

