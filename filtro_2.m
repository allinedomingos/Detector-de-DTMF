%% AE2 - Projeto FIR
% Grupo I: Alline Domingos, Fabiano Kraemer, Natália Miranda

%% Filtro II - Janela ajustável Kaiser(HP)
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%especificações
fa = 8000; % kHz 
Gp = 0;    % dB,
As = 30;   % dB
Ap = 0.5;  % dB  
fp = 950;  % Hz; 
fs = 1300; % Hz, 

%normalização:

wp = (2 * pi * fp) / fa; %Omega de passagem, a frequência de passagem do filtro HP
ws = (2 * pi * fs) / fa; %Omega de corte, a frequência de corte
%Lambdas
lp = 2 * fa * tan((wp/2));
ls = 2 * fa * tan((ws/2));
%Omegas
Os = lp / ls;
Op = wp / wp;


fcuts = [fp fs];
mags = [0 1];
%Os valores abaixo, Dstop e Dpass foram obtidos pelo FDATool
Dstop = 0.031622776602;  % Stopband Attenuation
Dpass = 0.028774368332;  % Passband Ripple
devs = [Dstop Dpass];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% para α, β e Δω, estimar esses valores utilizando a função kaiserord:
[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fa)
%n = 40; testes

%a função fir1 usa o método janela para gerar o filtro, também foi
%adicionado um pequeno valor ao beta para o equiripple não passar do 0.25
h_fir = fir1(n,Wn,'high',kaiser(n+1,beta+0.04),'noscale');

den = 1; %denominador
num = real(h_fir);

%% Plotando resposta em freq
[Hw,w] = freqz(h_fir);
hold on;
grid on;
figure(1);
plot(w*fa/2/pi,20*log10(abs(Hw)))
title(['Kaiser filter N = ' num2str(n)])
ylim([-As-10 2]);
xlim([500 2500]);
% mascaras
plot([0 fp fp fa/2],[-As -As Ap/2 Ap/2], ':r');
plot([fs fs fa/2],[-As -Ap/2 -Ap/2], ':m');

title(['n = ', num2str(n) ', Resposta em Frequência'] );

hold off;

%% Plotando polos e zeros

% Polos e zeros

%visualizando o plano z de H(p)
hold on;
grid on;
figure(3);
zplane(num,den);
title('Polos e zeros');
hold off;


%% Resposta em fase

grid on;
hold on;
figure(4);
phasez(h_fir);
ylim([-As-10 10]);


%% 9 PASSO - Atraso de Grupo
figure(4)
grpdelay(h_fir);
title('Atraso de grupo de H(z)');

%%
fvtool(h_fir,1) % usado para conferência dos gráficos