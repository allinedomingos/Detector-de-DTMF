%% AE2 - Projeto FIR
% Grupo I: Alline Domingos, Fabiano Kraemer, Natália Miranda

%% Filtro I - Butterworth (LP)
clear all;
clc;
close all;

% Especificações
fp = 950;  % freq. de passagem
fs = 1300; % freq. de rejeição
fa = 8000; % freq. de amostragem
As = 30;   % Atenuação na freq. rejeição
Ap = 0.5;  % Atenuação na freq. de passagem
Gp = 0;    % Ganho em db na banda de passagem

%% 1 PASSO - Normalização das frequencias
Wp = fp/(fa/2);
Ws = fs/(fa/2);
%d_w = -0.035;
d_w = 0;
Wc = ((Wp+Ws)/2 + d_w)*pi;
%d_w1 = abs(Wp - Ws)*pi;

%% 2 PASSO  -  Estimativa da ordem do filtro
M = 32;

%% 3 PASSO - Verificar M: onde M é N / 2 para N par e (N + 1) / 2 para N impar
n = -M:M;
%% 4 PASSO - Calculo dos coeficientes da série de fourier
c_hp = (sin(Wc*n)./(pi*n));
c_hp(M+1) = (Wc./pi);

%% 5 PASSO - Calculo janela
% Bartlett-Hanning  ou Hann
%w = barthannwin(2*M+1);
w = hann(2*M+1);
%w = bartlett(2*M+1);
%w = hamming(2*M+1);
%w =  blackman(2*M+1);
%resposta ao impulso
w = w';
b = c_hp.*w;
%% 6 PASSO - Polos e zeros
figure(1)
zplane(b,1);
title('Pólos e zeros');
[h,w]=freqz(b,1);
%title(['n = ' num2str(length(c_hp)-1)])

%% 7 PASSO - Filtro LP
figure(2)
plot((w/pi)*(fa/2), 20*log10(abs(h)))
hold on
grid on
%plot([fp fs],[-Ap -As],'k*');
plot([0 fs fs fs],[Ap/2 Ap/2 -As -As],'--m');
plot([0 fp fp],[-Ap/2 -Ap/2 -As],'--m');
ylim([-As 2]);
hold off
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
title('Filtro H(z)');

%% 8 PASSO - Resposta em frequência 
figure(3)
freqz(b,1);
title('Resposta em Frequencia de H(z)');

%% 9 PASSO - Atraso de Grupo
figure(4)
grpdelay(b);
title('Atraso de grupo de H(z)');
