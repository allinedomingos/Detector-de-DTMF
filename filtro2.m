%% AE1 - Projeto IIR
% Grupo I: Alline Domingos, Fabiano Kraemer, Natália Miranda
clear all;
close all;
clc;

%% Filtro II - Chebychev I (HP)
% Especificações
fp = 1300;  % Hz Frequência de passagem
fs = 950; % Hz Frequência de corte
As = 30;   % dB - atenuação em fs

%valores default
fa = 8000; % Hz Numero de amostras
Gp = 0;	% dB Ganho na banda de passagem
GdB = 10^(Gp / 20);
Ap = 0.5;  % dB - atenuação na Fp


%% 1) Passo, normalização:

wp = (2 * pi * fp) / fa; %Omega de passagem, a frequência de passagem do filtro LP
ws = (2 * pi * fs) / fa; %Omega de corte, a frequência de corte

%Lambdas
lp = 2 * fa * tan((wp/2));
ls = 2 * fa * tan((ws/2));

%Omegas
Os = lp / ls;
Op = wp / wp;

%% 2) Passo, achando a ordem do filtro, polos e zeros:

% n para filtro Chebyshev:
n = (acosh(sqrt(((10^(0.1 * As)) - 1) / ((10^(0.1 * Ap)) - 1)))) / (acosh(Os));
n = ceil(n);

epson = sqrt((10^(0.1 * Ap)) -1);
r = epson^(-(1/n));

%Cálculo de PK
k = 1:n;
wk = ((2 * k - 1) * pi) / (2.* n);
fi2 = (1/n) * asinh(1/epson);
pk = -sinh(fi2) * sin(wk) + 1j*(cosh(fi2)*cos(wk)); % polos do filtro

% Ajuste do ganho
if mod(n,2) == 0
	Gp = sqrt(1/(1 + epson ^2));   
else Gp = 1;
end

% Polos e zeros
den = real(poly(pk)); %denominador
num = den(end);


%visualizando o plano z de H(p)
figure(1);
zplane(num,den);
title('Polos e zeros do Prototipo H(p)');


%% Cálculo do protótipo com funções simbólicas
num = num * 0.944;
[h,w] = freqs(num, den, 2000);
hdB = 20*log10(abs(h));

figure(2);
grid on;
hold on;
ylim([-32, 2]);
xlim([0 2]);
plot(w, hdB);
plot([Op, Os],[-Ap -As],'k*');
plot([0 Op Op], [-Ap -Ap -100], '--r');
plot([0 Os Os Os*1.5], [0 0 -As -As], '--r');
xlabel('Omega p normalizado');
ylabel('Magnitude (dB)');
title('Magnitude do Prototipo H(p)');
hold off;

%% II Etapa - Passado H(p) para H(s)

syms s;
syms p;
inversao = lp / s;

Np = poly2sym(num, p);
Dp = poly2sym(den, p);

Hp(p) = Np/Dp; %

Hs(s) = subs(Hp(p), inversao); %

[Ns, Dens] = numden(Hs(s)); %
bs = sym2poly(Ns);
as = sym2poly(Dens);

figure(3);
freqs(num, den, 2000);
title('Resposta em frequência do Protótipo H(p)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');

%% II Etapa - Passado H(p) para H(s)
syms s;
D_s1(s) = poly2sym(den/num,s); 
D_s(s)= D_s1(s/lp);

Hsc(s) = collect(Hs(s));


[Nums, Dens] = numden(Hsc(s));
bs = sym2poly(Nums);
as = sym2poly(Dens);

figure(4)
hold on; 
grid on;
[h, w] = freqs(bs, as, 2000);
plot(w, 20*log10(abs(h))); ylim([-As-20 10]); 
%plot([fp, fs],[-Ap -As],'k*');
plot([-As ls ls lp+2000], [-As -As 0 0], '--r');
plot([-As-Ap lp lp lp+2000], [-35 -35 -Ap -Ap], '--r');
title('Filtro passa faixa');
xlabel('Frequências do Filtro em Lambda (Hz)');
ylabel('Magnitude (dB)');
axis([0 lp+2000 -35 2]);
title('Magnitude do Filtro Analógico H(s)');
hold off;

figure(5);
zplane(bs,as);
title('Pólos e zeros do Filtro Analógico H(s)');


figure(6);
freqs(bs,as);
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
H_z(z) = collect(subs(Hsc,s));

% Polos e zeros
[num_z den_z] = numden(sym(H_z));
bz = sym2poly(num_z);
az = sym2poly(den_z);

w_hz = 0:0.01:pi;
calc_z = H_z(exp(1i*w_hz));

figure(7);
hold on;
plot([fp fs],[-Ap -As],'k*');
plot([-As fs fs 2000], [-As -As 0 0], '--r');
plot([fp fp fp 2000], [-As-Ap -As-Ap -Ap -Ap], '--r');
plot(fa*w_hz/pi/2,20*log10(abs(calc_z)));
xlabel('Frequências do Filtro (Hz)');
ylabel('Magnitude (dB)');
axis([0 1500 -32 2]);
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
