clc;
clear; 
close all;  

dt=9.629166007732353e-13;

load('XGD_EH_RIPML_Debye_ADE_1_77.mat')
load('XGD_EH_RIPML_Debye_RI_1_77.mat') 

figure(6)
plot(1e9*dt*(1:9347),9.5*XGD_EH_RIPML_Debye_ADE_1_77/max(abs(XGD_EH_RIPML_Debye_ADE_1_77)),'-b','color','[0.08 0.17 0.55]','linewidth',3)
hold on
plot(1e9*dt*(1:9347),9.5*XGD_EH_RIPML_Debye_RI_1_77/max(abs(XGD_EH_RIPML_Debye_RI_1_77)),'--r','linewidth',3)
axis([0 9 -8 10])
ylabel('Electric Field Ez (V/m)')
xlabel('Time (ns)')
legend('FDTD-RI method','FDTD-ADE method')


load('RI_SAR_1')
load('ADE_SAR_1')

figure(7)
imagesc(Brain_SAR_single_RI_3_2D(end:-1:1,:)/max((max(Brain_SAR_single_RI_3_2D))));
colormap(jet);
colorbar
axis([0 300 0 300])
ylabel('(mm)')
xlabel('(mm)')


figure(77)
imagesc(Brain_SAR_single_ADE_3_2D(end:-1:1,:)/max((max(Brain_SAR_single_ADE_3_2D))));
colormap(jet);
colorbar
axis([0 300 0 300])
ylabel('(mm)')
xlabel('(mm)')


% Grid for EM fields
Imax = 298;  
Jmax = 298; 
Kmax = 170; 

% PML thickness in each direction 
NPML = 11; 

fp = fopen('Head_voxel.dat');
data=fread(fp,'int8');
Brain_data=reshape(data,256,256,128);
Brain_p=zeros(Imax+1, Jmax+1, Kmax+1);

for i = 1:256
	for j = 1:256
	    for k = 1:128
		    Brain_p(i+NPML+5,j+NPML+5,k+NPML+10) =  Brain_data(i,j,k);
        end
    end
end

%% RI SAR
figure(8)
isosurface(Brain_p(:,:,end:-1:1),0);
alpha(0.2)
hold on;
isosurface(Brain_p(:,:,end:-1:1),110,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),123,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),121,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),120,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),119,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),118,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),117,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),116,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),115,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),114,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),113,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),112,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),111,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),109,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),108,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on
isosurface(Brain_p(:,:,end:-1:1),107,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),106,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),105,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),104,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on
isosurface(Brain_p(:,:,end:-1:1),103,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),102,300*Brain_SAR_total_RI_4/max(max(max(abs(Brain_SAR_total_RI_4)))));
ylabel('y')
xlabel('x')
zlabel('z')
colormap(jet);
colorbar
caxis([0,1])

%% ADE SAR
figure(88)
isosurface(Brain_p(:,:,end:-1:1),0);
alpha(0.2)
hold on;
isosurface(Brain_p(:,:,end:-1:1),110,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),123,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),121,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),120,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),119,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),118,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),117,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),116,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),115,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),114,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),113,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),112,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),111,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),109,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),108,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on
isosurface(Brain_p(:,:,end:-1:1),107,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),106,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),105,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),104,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on
isosurface(Brain_p(:,:,end:-1:1),103,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
hold on;
isosurface(Brain_p(:,:,end:-1:1),102,282*Brain_SAR_total_ADE_4/max(max(max(abs(Brain_SAR_total_ADE_4)))));
ylabel('y')
xlabel('x')
zlabel('z')
colormap(jet);
colorbar
caxis([0,1])












