%% AE1 - Projeto IIR
% Grupo I: Alline Domingos, Fabiano Kraemer, Natalia Miranda
% Exemplo base aula 15 10/04/18
clear all; close all;clc;
%% Filtro III - Eliptico (BP1)
% Especificacações
fa = 8000; % freq. de amostragem

fs1 = 627; % freq. de rejeição
fp2 = 683; % freq. de passagem
fp3 = 711; % freq. de passagem
fs4 = 767; % freq. de rejeição
G0 = 0;% Ganho em db na banda de passagem
Ap = 0.5; %Atenuação na freq. de passagem 
As = 30; % Atenuação na freq. rejeição   

%% I Etapa - projeto de um filtro passa faixa (BP) protótipo normalizado H(p) com frequência de passagem Ωs = 1
% Normalizando as frequências
% Frequencia angular 
ws1 = (2*pi*fs1/fa); % divide por fa para evitar a distorção
wp2 = (2*pi*fp2/fa);
wp3 = (2*pi*fp3/fa);
ws4 = (2*pi*fs4/fa);
w0 = sqrt(fp2*fp3)*2*pi;

% Lambda
ls1 = 2*fa*tan(ws1/2);
lp2 = 2*fa*tan(wp2/2);
lp3 = 2*fa*tan(wp3/2);
ls4 = 2*fa*tan(ws4/2);
l0 = sqrt(lp2*lp3);

Bl = lp3-lp2;
% Omega
Ws = abs((l0^2 - ls4^2)/(Bl*ls4));
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
f_s = (s^2 + l0^2)/(Bl*s); %p = (1/Bl)*((Ws^2) + (l0^2)/ Ws); 

pretty(vpa(f_s, 4));
Hp(p) = Np/Dp;
pretty(vpa(Hp(p), 4));
w_hp = linspace(0, Ws+1, fa);
Hp_p = 20*log10(abs(Hp(1j*w_hp)));

%Plotagem
figure(1)
plot(w_hp, Hp_p, 'b');
grid on;hold on;
plot([0 Wp Wp], [-Ap -Ap -100], '--r');
plot([0 Ws-1.9 Ws-1.9 Ws*1.5], [0 0 -As -As], '--r');
axis([0 Ws*1.2 -50 5 ]);
xlabel('Frequencia normalizada (\Omega)','fontsize',13); ylabel('Magnitude (dB)','fontsize',13);
title('H(p)','fontsize',14);

figure(2)
%zplane(b, a);
[nup, dep] = numden(Hp);
zplane(sym2poly(nup),sym2poly(dep));
title('Polos e Zeros de H(p)','fontsize',14);
xlabel('Real','fontsize',13); ylabel('Imagin?rio','fontsize',13);

% Fase
figure(3)
freqs(sym2poly(nup),sym2poly(dep))
title('Magnitude e Fase H(p)');

%% II Etapa - Passado H(p) para H(s)
Hs(s) = vpa(collect(Hp((s^2 + l0^2)/(Bl*s))));
pretty(vpa(Hs(s), 4));

w_hs = linspace(0.01, 5000*2*pi, fa);
Hs_s = 20*log10(abs(Hs(1j*w_hs)));

%Plotagem
figure(4)
plot(w_hs, Hs_s, 'b');
grid on;hold on;
plot([0 ls1 ls1 ls4 ls4 ls4*2], [-As -As 0 0 -As -As ],'--r');
plot([lp2 lp2 lp3 lp3 ], [-90 -Ap -Ap -90], '--r');
axis([ls1*0.5 ls4*1.2 -50 5 ])
xlabel('Frequencia (\lambda)','fontsize',13); ylabel('Magnitude (dB)','fontsize',13);
title('H(s)','fontsize',14);
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
plot([0 fs1 fs1 fs4 fs4 fs4*2], [-As -As 0 0 -As -As ],'--r');
plot([fp2 fp2 fp3 fp3 ], [-90 -Ap -Ap -90], '--r');
xlabel('Frequencia (Hz)','fontsize',13); ylabel('Magnitude (dB)','fontsize',13);
title(['Filter Bandpass - Ordem: ' num2str(ordem)],'fontsize',14);
axis([fs1-200 fs4+200 -50 5 ])

%zplane
figure(8)
[nuz, dez] = numden(Hz);
zplane(sym2poly(nuz),sym2poly(dez));
title('Polos e Zeros  de H(z)','fontsize',14);
xlabel('Real','fontsize',13); ylabel('Imagin?rio','fontsize',13);

%Fase
figure(9)
freqz(sym2poly(nuz),sym2poly(dez))
title('Magnitude e Fase H(z)');'H(z)'

%%  Atraso de grupo 
figure(10);
grpdelay(sym2poly(nuz),sym2poly(dez));
title('Atraso de grupo');