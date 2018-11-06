%% AE2 - Projeto FIR
% Grupo I: Alline Domingos, Fabiano Kraemer, Natalia Miranda
%O filtro BP deve ser projetado usando o algoritmo de Parks-McCleallan.
clear all; close all;clc;
%% Filtro III - Band Pass
% Especificacações
fa = 8000; % freq. de amostragem, fsamp

fs1 = 627; % freq. de rejeição
fp2 = 683; % freq. de passagem
fp3 = 711; % freq. de passagem
fs4 = 767; % freq. de rejeição
G0 = 0;% Ganho em db na banda de passagem
Ap = 0.5; %Atenuação na freq. de passagem 
As = 30; % Atenuação na freq. rejeição  

% Vetor de bordas de banda de freqüência
f = [fs1 fp2 fp3 fs4];
% Vetor contendo as amplitudes desejadas nos pontos especificados em f.
a = [0 1 0];
% Vetor do mesmo tamanho que especifica o desvio ou ondulações máximo 
%permitido entre a resposta de freqüência e a amplitude desejada do filtro 
%de saída para cada banda.
Dstop1 = 0.02;  % First Stopband Attenuation
Dpass  = 0.018;  % Passband Ripple
Dstop2 = 0.02;  % Second Stopband Attenuation
dev = [Dstop1 Dpass Dstop2]; 

%% Cálculo da ordem dos parâmetros usando FIRPMORD.
[n, fo, ao, wpm] = firpmord(f, a, dev, fa);
%% Cálculo do coeficiente usando a função FIRPM.
b = firpm(n, fo, ao, wpm, 'low');

%% Plotagem dos Gráficos
figure(1)
freqz(b,1,8000,fa)
title('Hz');
[hz,wz] = freqz(b,1,8000,fa);

figure(2)
plot(wz, 20*log10(abs(hz))); hold on; grid on;
title(['Filter de Ordem= ' num2str(n)]);
plot([0 fs1 fs1 fs4 fs4 fa/2], [-As -As Ap/2 Ap/2 -As -As], ':r');
plot([fp2 fp2 fp2 fp3 fp3 fa/2], [-As*100 -As*100 -Ap/2 -Ap/2 -As*100 -As*100], ':k' );
axis([500 900 -50 1 ])
xlabel('Frequencia (Hz)','fontsize',13); ylabel('Magnitude (dB)','fontsize',13);

figure(3);zplane(b,1);title('Pólos e zeros');
xlabel('Real','fontsize',13); ylabel('Imaginario','fontsize',13);

%%  Atraso de grupo 
figure(4)
grpdelay(b)
title('Atraso de grupo');