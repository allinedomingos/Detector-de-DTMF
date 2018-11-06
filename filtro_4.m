%% AE1 - Projeto IIR
% Grupo I: Alline Domingos, Fabiano Kraemer, Natália Miranda
clear all; close all;clc;

%% Filtro IV - Elíptico (BS1)
% Observação: Para projetar um filtro bandstop é preciso projetar um
% passa-baixa e transformar em bandstop. Passos seguidos do livro, páginas
% 216 à 218.
 
% Especificacações
fa = 8000; % freq. de amostragem

fp1 = 1202; % freq. de passagem
fs2 = 1309; % freq. de rejeição
fs3 = 1363; % freq. de rejeição
fp4 = 1470; % freq. de passagem

G0 = 0;% Ganho em db na banda de passagem
Ap = 0.5; % Atenuação na freq. de passagem 
As = 30; % Atenuação na freq. rejeição  
%% I Etapa - projeto de um filtro rejeita faixa (BS) protótipo normalizado H(p) com frequência de passagem Ωs = 1
% Normalizando as frequências
% Frequencia angular 
wp1 = (2*pi*fp1/fa); 
ws2 = (2*pi*fs2/fa);
ws3 = (2*pi*fs3/fa);
wp4 = (2*pi*fp4/fa);
% Lambda
lp1 = 2*fa*tan(wp1/2);
ls2 = 2*fa*tan(ws2/2);
ls3 = 2*fa*tan(ws3/2);
lp4 = 2*fa*tan(wp4/2);  

% Valor teorico - freq de passagem
l0 = sqrt(lp1*lp4);
w0 = sqrt(fp1*fp4)*2*pi;
Bl = lp4-lp1;

% Omega
Ws1 = abs((Bl*ls2)/(l0^2 - ls2^2));
Ws2 = abs((Bl*ls3)/(l0^2 - ls3^2));
Ws = min(Ws1, Ws2);
Wp = 1;

%% Eliptico
[n, Wn] = ellipord(Wp, Ws, Ap, As,'s');
[b,a] = ellip(n,Ap,As, Wn, 's');
ordem = 2*n;
%% Calculo do prototipo com funções simbolicas
syms p s z;

Np = poly2sym(b, p);
Dp = poly2sym(a, p);

% Função de transferencia em p
f_s = Bl*(s/(s^2 + w0^2)); 
pretty(vpa(f_s, 4));

Hp(p) = Np/Dp;
pretty(vpa(Hp(p), 4));
w_hp = linspace(0, Ws+1, fa);
Hp_p = 20*log10(abs(Hp(1j*w_hp)));

%Plotagem
figure(1)
plot(w_hp, Hp_p, 'b');
grid on;hold on;
plot([0 Wp Wp], [-Ap -Ap -107], '--r');
plot([0 Ws Ws Ws*2], [0 0 -As -As], '--r');
axis([0 6 -107 0])

xlabel('Frequencia normalizada (\Omega)','fontsize',13); ylabel('Magnitude (dB)','fontsize',13);
title('H(p)','fontsize',14);

%zplane
figure(2)
[nup, dep] = numden(Hp);
zplane(sym2poly(nup),sym2poly(dep));
title('Polos e Zeros de H(p)','fontsize',14);
xlabel('Real','fontsize',13); ylabel('Imagin?rio','fontsize',13);

% Fase
figure(3)
freqs(sym2poly(nup),sym2poly(dep))
title('Magnitude e Fase H(p)');
%% II Etapa - Passado H(p) para H(s)
Hs(s) = vpa(collect(Hp(Bl*(s/(s^2 +l0^2)))));
pretty(vpa(Hs(s), 4));

w_hs = linspace(0.01, 5000*2*pi, fa);
Hs_s = 20*log10(abs(Hs(1j*w_hs)));

%Plotagem
figure(4)
plot(w_hs, Hs_s, 'b');
grid on;hold on;
plot([0 lp1 lp1], -[Ap Ap As], '--r')
plot([0 ls2 ls2 ls3 ls3 lp4+100], [0 0 -As -As 0 0], '--m')
plot([lp4+1 lp4 lp4], -[Ap Ap As], '--r')
axis([ls2*0.5 ls3*1.2 -50 5 ])
xlabel('Frequencia (\lambda)','fontsize',13); ylabel('Magnitude (dB)','fontsize',13);
title('H(s)','fontsize',14);

%zplane
figure(5)
[nus, des] = numden(Hs);
zplane(sym2poly(nus),sym2poly(des));
title('Polos e Zeros de H(s)','fontsize',14);
xlabel('Real','fontsize',13); ylabel('Imagin?rio','fontsize',13);

%Fase
figure(6)
freqs(sym2poly(nus),sym2poly(des))
title('Magnitude e Fase H(s)');

%% III Etapa - Passando Hs(s) para Hz(z)
Hz(z) = vpa(collect(Hs(2*fa*(z-1)/(z+1)))); 
pretty(vpa(Hz(z), 4));
Hzc(z) = collect(Hz(z));
pretty(vpa(Hzc(z), 4));

w_hz = linspace(0.01, pi,fa);
Hz_z = 20*log10(abs(Hz(exp(1i*w_hz))));

%PLotagem
figure(7)
plot((w_hz)*(fa/(2*pi)), Hz_z, 'b');
grid on;hold on;
plot([fp1-100 fp1 fp1], -[Ap Ap As], '--r')
plot([fp1-100 fs2 fs2 fs3 fs3 fp4+100], [0 0 -As -As 0 0], '--m')
plot([fp4+100 fp4 fp4], -[Ap Ap As], '--r')
xlabel('Frequencia (Hz)','fontsize',13); ylabel('Magnitude (dB)','fontsize',13);
title(['Filter Bandpass - Ordem: ' num2str(ordem)],'fontsize',14);

%zplane
figure(8)
[nuz, dez] = numden(Hz);
zplane(sym2poly(nuz),sym2poly(dez));
title('Polos e Zeros  de H(z)','fontsize',14);
xlabel('Real','fontsize',13); ylabel('Imagin?rio','fontsize',13);

%Fase
figure(9)
freqz(sym2poly(nuz),sym2poly(dez))
title('Magnitude e Fase H(z)');

%%  Atraso de grupo 
figure(10);
grpdelay(sym2poly(nuz),sym2poly(dez));
title('Atraso de grupo');