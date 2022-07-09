clc;
close;
close all;


sam=2000;
Fre=linspace(0,5e9,sam);


load('ADE_T_fft_1');
load('RI_T_fft_1');
load('T_analytical');

figure(3)
plot(Fre(1:60:2000),10*log10(T_analytical(1:60:2000)),'ok')
hold on
plot(Fre(1:10:2000),10*log10(RI_T_fft_1(1:10:2000)),'-b','color','[0.08 0.17 0.55]','linewidth',3)
hold on
plot(Fre,10*log10(ADE_T_fft_1),'--r','linewidth',3)
axis([0 5e9 -1.4 0])
legend('Analytical Method','FDTD-RI Method','FDTD-ADE Method')
ylabel('Transmission Coefficient (dB)')
xlabel('Frequency (GHz)')


