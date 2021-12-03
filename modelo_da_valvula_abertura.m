%%
%AutoCERES
%Modelagem da válvula solenoide 
%20/08/2021 - Mariana Obiedo Piñeiro
%--------------------------------------
close all;
clear all;

%constantes do sistema
%subsistema eletromagnetico
R = 25.4; %resistencia do sistema RL
N = 1105; %numero de voltas do solenoide
de = 10e-3; %diametro do entreferro
Ae = pi*(de/2)^2; %area da seção do entreferro
Me = pi*4e-7; %permeabilidade magnetica do ar
Xe = 8e-3;%distancia em que desloca o embolo 8mm


%subsistema mecanico
%¨determinando ks, a constante elástica da mola, a partir do ensaio
xk = [2e-3 4e-3 6e-3 8e-3]; %deslocamento da mola em m
mk = [4e-3 8e-3 12e-3 16e-3]; %massa aplicaga em kg
g = 9.807; %gravidade 
Fk = g*mk;

%figure();
%plot(xk,Fk);
%grid();
%ks=polyfit(xk,Fk,1);
ks = (Fk/xk)

ms = 4e-3; %massa do embolo 4g

wn = sqrt(ks/ms);
ta = 14.67e-3;%14,67e-3;
zs = 3/(ta*wn); %ta é o tempo de acomodação
bs = 2*zs*wn*ms; %coeficiente de atrito viscoso

zeta = 12e-3; %pré tensão da mola 

%subsistema hidraulico
Cc = 0.66; %coeficiente de contração típico
Cv = 0.98; %coeficiente de velocidade típico
Cd = Cc*Cv; %coeficiente de descarga
rho = 1e3; %densidade agua
A1 = 85.75e-6; %area de entrada do fluido 85,75mm² em m²
A2 = 17.3e-6; %area de saida do fluido 17,3mm² em m²
p1=4e5; %pressao interna 4 bar
p2=1e5;

Fpres = (A1-A2)*p1; 

% Tempo de amostragem
ts = 1/100000;

% Condições iniciais
I(2) = 0; % Corrente
V(1) = 0; % tensão
Vx(2) = 0;% velocidade
M(1) = 0; % vazão
X(1) = 0; %posição
X(2) = 0; %posição
L(1) = 0;
L(2) = 0;

for k = 2:100000
    V(k) = 12; %degrau de tensão 12Vcc
    tempo(k) = k*ts; %conversão do tempo de amostragem em segundos
    
    % Equações da discretização da válvula
    L(k) = (N*N*Ae*Me)/(Xe - X(k));
    L(k+1) = ts*(L(k)*L(k)*Vx(k)/(N*N*Ae*Me)) + L(k);
    I(k+1) = I(k) + ts*(V(k)/L(k))  - I(k)*((R*ts)/L(k)) - I(k)*ts*((L(k)*Vx(k))/(N*Ae*Me));
    %Fm(k) = N*I(k);
    Fm(k) = 0.5*I(k)*I(k)*((L(k+1)-L(k))/ts);
    %Fm(k) = (I(k)*I(k)*L(k)*L(k))/(N*N*Ae*Me);
    Vx(k+1) =  Vx(k) - ((ts*bs)/ms)*Vx(k) - ((ts*ks)/ms)*(X(k) - zeta) + (ts*Fm(k))+ (ts*Fpres);
    %Vx(k+1) =  Vx(k) - ((ts*bs)/ms)*Vx(k) - ((ts*ks)/ms)*(X(k) + zeta) + (ts*Fm(k))/ms + (ts*Fpres)/ms;
    
    if X(k) == 0
        M(k+1) = 0;
    else
        M(k+1) = ts*rho*Cd*A1*sqrt(((2*(p1-p2))/(rho*((A1/A2)-1)))) + M(k);
    end
    
    %a posição x(t) derivada da velocidade:
    X(k+1) = X(k) + Vx(k)*ts;
    
    if X(k+1) > Xe
        X(k+1) = Xe;
        X(k) = Xe;
    end
    if X(k+1) < 0
        X(k+1) = 0;
        X(k) = 0;
    end
    if I(k+1) > .432
            I(k+1) = .432;
            I(k) = .432;
    end

end

%transformação do tempo discreto para tempo contínuo
tempo(k+1) = (k+1)*ts;
I(k+1) = I(k);
X(k+1) = X(k);
L(k+1) = L(k);
Vx(k+1) = Vx(k);

figure;
plot(tempo,M,'r');
grid;
title ('X(t)') ;
xlabel('t [s]');
ylabel('x [m]');

