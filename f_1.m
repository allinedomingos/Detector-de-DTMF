%% AE1 - Projeto IIR
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

%% I Etapa - projeto de um filtro passa baixas (LP) protótipo normalizado H(p) com frequência de passagem Ωs = 1
% Normalizando as frequências
% Frequencia angular 
Wp = (2*pi*fp)/fa;
Ws = (2*pi*fs)/fa;
% Lambda 
lp = 2*fa*tan((Wp/2)); 
ls = 2*fa*tan((Ws/2));
% Omega
ohm_p = 1;
ohm_s = ls/lp; 
 

%% Determinando ordem, polos e zeros
% Ordem do filtro
n = log10(((10^(0.1*As)) - 1)/((10^(0.1*Ap)) -1))/(2*log10(ohm_s)); 
n = ceil(n)
epson = sqrt((10^(0.1*Ap)) -1);
r = epson^(-(1/n));

% Calculo de pk
k = 1:n;
p_k = r*exp(j*pi*((2*k+(n-1))/(2*n))); 

%Polos e zeros
den = real(poly(p_k));
num = den(end);

%Calculo simbolico
syms p ;
Den_p(p) = poly2sym(den/num,p);

% Calculo de H_0 - Example 4.3 pag 201
H_0 = Den_p(j*0);
num_p = H_0*(10^(Gp/20));
zeros_p = num_p

%% Calculo do prototipo com funções simbolicas
% Função de transferencia em p
H_p(p) = symfun(num_p/Den_p,p);
funcao_hp = H_p(p);

% Zplane de H(p)
figure(1);
zplane(num,den);
title('Pólos e zeros do Protótipo H(p)');

w_hp=0:0.01:ohm_s;
calc_p = H_p(j*w_hp); %jw

figure(2);
hold on;
plot(w_hp,20*log10(abs(calc_p)));
plot([ohm_p ohm_s],[-Ap -As],'k*');
plot([0 ohm_p ohm_p ohm_p],[-Ap -Ap -As -As*2],'--m');
plot([0 ohm_s ohm_s ohm_s],[0 0 -As -As*2],'--m');
hold off;
axis([0 1.5 -30 5]);
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
title('Magnitude do Protótipo H(p)');

np_aux = sym2poly(H_0);
dp_aux = sym2poly(Den_p(p));

figure(3);
freqs(np_aux,dp_aux)
title('Resposta em frequência do Protótipo H(p)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');

%% II Etapa - Passado H(p) para H(s)
syms s;
H_s(s) = H_p(s/lp); % calculo Hs
D_s1(s) = poly2sym(den/num,s); 
D_s(s)= D_s1(s/lp);

%w_hs=0:1e3:ls+lp; 
w_hs = logspace(0,log10(ls),1000); 
calc_s = H_s(1i*w_hs); %jw 
funcao_hs = H_s(s);

% Polos e zeros H(s)
den_s = sym2poly(D_s);
num_s = sym2poly(H_0);

figure(4);
hold on;
plot([lp ls],[-Ap -As],'k*');
plot([lp lp 0],[-Ap*100 -Ap -Ap],'--m');
plot([0 ls ls ls*10],[0 0 -As -As],'--m'); 
plot(w_hs,20*log10(abs(calc_s)));
xlabel('Frequências do Filtro em Lambda (Hz)');
ylabel('Magnitude (dB)');
axis([0 10000 -30 5]);
title('Magnitude do Filtro Analógico H(s)');
hold off;

figure(5);
zplane(num_s,den_s);
title('Pólos e zeros do Filtro Analógico H(s)');


figure(6);
freqs(num_s,den_s);
title('Resposta em frequência do Filtro Analógico H(s)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');


%% III Etapa - Passando Hs(s) para Hz(z)
%Método de tranformação bilinear
%A região à esquerda do eixo imaginário no Plano S é mapeada dentro do círculo de raio unitário no Plano Z.
%Calculo simbolico
syms z;
s = (2*fa*(z-1)/(z+1));
H_z(z) = collect(subs(H_s,s));
funcao_hz = vpa(H_z(z));

% Polos e zeros
[num_z den_z] = numden(sym(H_z));
bz = sym2poly(num_z);
az = sym2poly(den_z);

w_hz = 0:0.01:pi;
calc_z = H_z(exp(1i*w_hz));

figure(7);
hold on;
plot([fp fs],[-Ap -As],'k*');
plot([fp fp 0],[-Ap*100 -Ap -Ap],'--m'); 
plot([0 fs fs 10*fs],[0 0 -As -As],'--m'); 
plot(fa*w_hz/pi/2,20*log10(abs(calc_z)));
xlabel('Frequências do Filtro (Hz)');
ylabel('Magnitude (dB)');
axis([500 1300 -30 5]);
title('Magnitude do Filtro Digital H(z)');
hold off;

figure(8);
zplane(bz,az);
title('Pólos e zeros do Filtro Digital H(z)');

figure(9);
freqz(bz,az);
xlim([0 0.9]);
title('Resposta em frequência do Filtro Digital H(z)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
xlim([0 0.7]);
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');

%% Atraso de grupo
figure(10);
grpdelay(bz,az)
title('Atraso de grupo H(z)');