clear all; close all; clc;
%%___________________________________________________
%                  CARACTERIZAÇÃO LASER JDSU        |
%                          980nm                    |
%___________________________________________________|
%% DEFINIÇÃO DE PARÂMETROS

load('data_laser.mat');

%___________________________________________________%

%% CÁLCULOS

%Potência em mW
measured_data(:,3)=10.^((measured_data(:,2)+att)/10);

%% Figura 1 - Espectro

figure;
plot(osa_data(:,1),osa_data(:,2),'LineWidth',2)
%legend('Espectro','Location','SouthEast')
ylabel('Potência de Saída (dBm)');
xlabel('Comprimento de onda (nm)');
title('Espectro Laser JDSU 980nm - I_b_o_m_b_e_i_o = 60mA - Data:21/10/2015')
axis([974 980 0 1.2]);
grid on;

%saveas(gca,' name ','epsc') %gera versão renderizada do plot
%saveas(gcf,'Figuras/Figura1.bmp','bmp')

%% Figura 2 - Corrente x Potência (mW)

figure;
plot(measured_data(:,1),measured_data(:,3),'o');
hold on;
plot(measured_data(:,1),measured_data(:,3),'k','LineWidth',2);
ylabel('Potência de Saída (mW)');
xlabel('Corrente de Bombeio (mA)');
title('Caracterização Laser JDSU 980nm - Data:21/10/2015')
set(gca,'XTick',0:10:150,'YTick',0:5:80);
%axis([974 980 0 1.2]);
grid on;

%saveas(gcf,'Figuras/Figura2.bmp','bmp')

%% Finaliza (Limpa variáveis e Fecha figuras)
clear all
close all