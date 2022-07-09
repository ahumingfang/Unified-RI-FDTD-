clc;
clear;
close all;

load('a5')
load('a6')
load('b5')
load('b6')
load('A_RI_1')
load('A_RI_2')
load('A_ZT_1')
load('A_ZT_2')
load('A_ADE_1')
load('A_ADE_2')


dt = 0.34e-12;
Ntimesteps =3530; 



figure(3)
plot(a5,b5/max(abs(b5)),'ok')
hold on
plot(dt*(1:Ntimesteps)*1e9,A_RI_1(1:Ntimesteps)/max(abs(A_RI_1(1:Ntimesteps))),'-b','color','[0.08 0.17 0.55]','linewidth',3)
hold on
plot(dt*(1:Ntimesteps)*1e9,A_ZT_1(1:Ntimesteps)/max(abs(A_ZT_1(1:Ntimesteps))),'-.g','linewidth',3)
hold on
plot(dt*(1:Ntimesteps)*1e9,A_ADE_1(1:Ntimesteps)/max(abs(A_ADE_1(1:Ntimesteps))),'--r','linewidth',3)
set(gca,'xtick',0:0.2:1.2);
axis([0 1.2 -0.8 1])
legend('Reference [11]','FDTD-RI method','FDTD-ZT method','FDTD-ADE Method')
ylabel('Electric Field Ez (V/m)')
xlabel('Time (ns)')






figure(4)
plot(a6,b6/max(abs(b6)),'ok')
hold on
plot(dt*(1:Ntimesteps)*1e9,A_RI_2(1:Ntimesteps)/max(abs(A_RI_2(1:Ntimesteps))),'linewidth',3)
hold on
plot(dt*(1:Ntimesteps)*1e9,A_ZT_2(1:Ntimesteps)/max(abs(A_ZT_2(1:Ntimesteps))),'-.g','linewidth',3)
hold on
plot(dt*(1:Ntimesteps)*1e9,A_ADE_2(1:Ntimesteps)/max(abs(A_ADE_2(1:Ntimesteps))),'--r','linewidth',3)
set(gca,'xtick',0.2:0.2:1.2);
axis([0.2 1.2 -1 0.8])
legend('Reference [11]','FDTD-RI method','FDTD-ZT method','FDTD-ADE Method')
ylabel('Electric Field Ez (V/m)')
xlabel('Time (ns)')















