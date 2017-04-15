% Turbojet Ideal
% aqui corremos em pi_c

%zerando------------------------------------------------------------------%
clear variables
close all
clc
format long
%-------------------------------------------------------------------------%

% Inputs------------------------------------------------------------------%

m0    = 0:0.5:3;                 %M0 CONJUNTO FECHADO ()
M0    = m0';                     %transposta 
T0    = 216.7;                   %temperatura inicial    |K|
hpr   = 42798.4;                 %      |kJ/kg|
yc    = 1.4;                     %      |adimensional|
Cpc   = 1.004832;                %      |kJ/kgK|
yt    = 1.35;                    
Cpt   = 1.0969416;               %      |kJ/kgK|
yAB   = 
Tt4   = 1666.7;                  %temperaturas no fim da camara de combustao |K|
p0_p9 = 1;                       %      |adimensional|
pi_dmax = 0.98;                  %      |adimensional|
pi_c  = 2:2:40;                  %varredura em pi_c - entra nas colunas das matrizes
pi_b  = 0.98;                    %      |adimensional|
pi_n  = 0.99;                    %      |adimensional|
e_c   = 0.92;                    %      |adimensional|
e_t   = 0.91;                    %      |adimensional|
n_b   = 0.99;                    %      |adimensional|
n_m   = 0.98;                    %      |adimensional|


j = length(M0);%cada linha do grafico CONJUNTO FECHADO
n = length(pi_c);
%-------------------------------------------------------------------------%

%Outputs--------------------------------------------------------------------%
%       F/m0ponto = |N/(kg/s)|      -> empuxo
%       f         = |adimensional|  -> razao combustivel/ air
%       s         = |(mg/s)/N|      -> consumo espeficico de combustivel por empuxo
%       nT        = |adimensional|  -> eficiencia termica
%       nP        = |adimensional|  -> eficiencia propulsiva
%       n0        = |adimensional|  -> eficiencia total
%       nc        = |adimensional|  -> eficiencia
%       nt        = |adimensional|  -> eficiencia
%
%-------------------------------------------------------------------------%

%pre-alocando matrizes com zeros para preenche-las------------------------%
tau_r   = zeros(j,n);
tau_lambda = zeros(j,n);
raz_1   = zeros(j,n);
raz_2   = zeros(j,n);
f       = zeros(j,n);
s       = zeros(j,n);
nP      = zeros(j,n);
nT      = zeros(j,n);
n0      = zeros(j,n);
%-------------------------------------------------------------------------%



%loop de varredura em pi_c--------------------------------------------------%
for i = 1:j
    for k = 1:n

    Rc          = ((yc-1)/yc)*Cpc;                        %1
    Rt          = ((yt-1)/yt)*Cpt;                        %2
    a0          = sqrt(yc*Rc*1000*T0);                     %3
    V0          = a0*M0;                                  %4
    tau_r       = 1+((yc-1)/2).*((M0).^2);                %5
    pi_r        = (tau_r).^(yc/(yc-1));                   %6
        if M0 <= 1
              n_r = 1;
        else  n_r = 1 - 0.075.*((M0 - 1).^1.35);          %7
        end 
    pi_d        = pi_dmax.*n_r;                           %8
    tau_lambda = (Cpt*Tt4)/(Cpc*T0);                      %9
    tau_c       = (pi_c).^((yc-1)/(yc*e_c));              %10
    nc          = (pi_c.^((yc-1)/(yc))-1)/(tau_c-1);      %11
    f           = (tau_lambda-tau_r.*tau_c)./(((n_b*hpr)/(Cpc*T0))-tau_lambda); %12
    tau_t       = 1-(1./(n_m.*(1+f))).*(tau_r./tau_lambda).*(tau_c-1);         %13
    pi_t        = tau_t.^(yt./((yt-1).*e_t));              %14
    nt          = (1-tau_t)/(1-tau_t.^(1/e_t));           %15
    pt9_p9      = (p0_p9).*pi_r.*pi_d.*pi_c.*pi_b.*pi_t.*pi_n;                %16
    M9          = sqrt((2/(yt-1)).*((pt9_p9).^((yt-1)/(yt))-1));              %17
    T9_T0       = ((tau_lambda.*tau_t)./((pt9_p9).^((yt-1)/yt))).*(Cpc/Cpt);    %18
    V9_a0       = M9.*sqrt((yt./yc).*(Rt./Rc).*(T9_T0));   %19
    F_m0ponto   = a0.*((1+f).*V9_a0-M0+(1+f).*(Rt/Rc).*(T9_T0./V9_a0).*((1-p0_p9)./yc)); %20
    s           = (f./(F_m0ponto)).*1000000;                          %21
    nT          = ((a0.^2).*((1+f).*(V9_a0).^2-(M0.^2)))./(2.*f.*hpr.*1000); %2
    nP          = (2.*V0.*(F_m0ponto))./((a0.^2).*((1+f).*(V9_a0).^2-M0.^2)); %23
    n0          = nP.*nT;                                 %24

    %Graficos---------------------------------------------------------%
        figure(1)
        subplot(2,2,1)
        plot(pi_c,F_m0ponto)
        ylim([0 1200])
        xlabel('pi_c')
        ylabel('F/m_0')
        grid
        grid minor
        legend('0','0.5','1','1.5','2','2.5','3')

        subplot(2,2,2)
        plot(pi_c,f)
        ylim ([0 0.045])
        xlabel('pi_c')
        ylabel('f')
        grid
        grid minor

        subplot(2,2,3)
        plot(pi_c,s)
        ylim ([20 65])
        xlabel('pi_c')
        ylabel('s')
        grid
        grid minor

        subplot(2,2,4)
        plot(pi_c,nT,pi_c,nP,pi_c,n0)
        ylim([0 1])
        xlabel('pi_c')
        ylabel('n')
        grid
        grid minor
        legend('nT','nP','n0')
    end
end
       
%----------------------------------------------------------------------%
