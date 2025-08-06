
%%Arbeitspunkt des kontinuierlichen nichtlinearen systems
xr1 =  0.18;
xr2 =  0.15;
xr3 =  0.12;
xR = [xr1;xr2;xr3];
ur1 =  5*10^-5;
ur3 =  10^-4;
uR = [ur1;ur3];
T = 0.1;
vxr1 = (xr1-0.31795)/-0.3314;
vxr3 = (xr3-0.30988)/-0.33781;
vxr = [vxr1;vxr3;vxr1;vxr3];






%% Dreitank Darstellung

% 0. Parameter
para.At = 0.0154;
para.Av = pi*(0.0040)^2;
para.a12 = 0.3485;
para.a10 = 0.6677;
para.g = 9.80665;
para.a23 = 0.3485;
para.a30I = 0.6591;
para.Ts = 0.1;

syms a10 a12 a23 a30I At Av g u1 u2 Ts
x        = sym('x',[3,1]); % Zustandsvektor
y        = sym('y',[3,1]); % Ausgangsvektor
u        = sym('u',[2,1]); % Ausgangsvektor

dx = [(u(1) - a12 * Av * sign(x(1) - x(2)) * sqrt(2*g*abs(x(1) - x(2))) - a10 * Av * sqrt(2*g*x(1))) / At;
    (a12 * Av * sign(x(1) - x(2)) * sqrt(2*g*abs(x(1) - x(2))) - a23 * Av * sign(x(2)-x(3)) * sqrt(2*g*abs(x(2)-x(3)))) / At;
    (u(2) + a23 * Av * sign(x(2)-x(3)) * sqrt(2*g*abs(x(2)-x(3))) - a30I * Av * sqrt(2*g*x(3))) / At];



%% Diskrete Übergangsfunktion bestimmen mit vorwärts Euler
simulate_step = @(x, u, T, para) ...
    ode45_step(@(t, x_) dreitank_model(x_, u, para), x, T);

function x_next = ode45_step(f, x0, T)
    [~, x_out] = ode45(f, [0 T], x0);
    x_next = x_out(end,:)';
end




%% Neuer Arbeitspunkt: numerische Bestimmung von xR mit f_d(xR, uR) = xR/Gleichung lösen

% Funktion f_d wie ist diskretes nicht lineares Modell
f_d = @(x, u) simulate_step(x, u, T, para);  % simulate_step wie vorher definiert

% esn wird gesucht xR für f_fix(x) = f_d(x, uR) - x = 0
f_fix = @(x) f_d(x, uR) - x;

% Startwert xR des kontinuierlichen Systems
xR_guess = xR;

% fsolve zur Bestimmung des Fixpunkts
options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-10);
[xR_new, fval, exitflag] = fsolve(f_fix, xR_guess, options);

if exitflag <= 0
    warning('fsolve hat keinen gültigen Arbeitspunkt gefunden.');
else
    disp('Neuer Arbeitspunkt xR:');
    disp(xR_new);
end

% neuer Arbeitspunkt festlegen
xR = xR_new;
vxr = [xR(1);xR(3);xR(1);xR(3)];

 %% Numerische Linearisierung um Arbeitspunkt
 delta = 1e-6;
 n = length(xR); m = length(uR);
 A_d = zeros(n);
 B_d = zeros(n, m);
 
 % A-Matrix/Dynamikmatrix bestimmen
 for i = 1:n
     dx = zeros(n,1); dx(i) = delta;
     f1 = f_d(xR + dx, uR);
     f2 = f_d(xR - dx, uR);
     A_d(:,i) = (f1 - f2) / (2*delta);
 end
 
 % B-Matrix/Konstantematrix bestimmen
 for i = 1:m
     du = zeros(m,1); du(i) = delta;

     f1 = f_d(xR, uR + du);
     f2 = f_d(xR, uR - du);
     B_d(:,i) = (f1 - f2) / (2*delta);
 end




%% Berechnung der linearisierten Matrizen für das ESN
nutzung_esn;
if esn.use_bias && esn.use_bias_out == false
    disp("mit input bias");
    eingang1 = vxr .* ones(4,50000);
    um = [vxr;esn.bias*esn.input_scaling];%arbeitspunkt Tankstände
elseif esn.use_bias && esn.use_bias_out
    disp("mit beiden bias")
    eingang1 = vxr .* ones(4,50000);
    um = [vxr;esn.bias*esn.input_scaling];
else
    disp("ohne input bias");
    eingang1 = vxr .* ones(4,50000);
    um = vxr;
end
states = esn.run(eingang1);
statem = states(:,49900);
xm = statem;%arbeitspunkt interne zustände des ESN

Wres = esn.Wres;
Win = esn.Win;

alpha = 0.75;


[C, D] = linearizeESN(Wres, Win, xm, um, alpha);


%%Berechnung geschlossener Regelkreis
Wout = esn.Wout(1:2,1:100);
    H = [1 0 0 ; 0 0 1 ]; 
    A_tilde = AN;%Matrix des linearen Systems
    B_tilde = bN; %matrix des linearen Systems
    k1 = [1/-0.3314 0 ; 0 1/-0.33781 ];% für umrechnung der Füllstände in volt bereich
    k2 = [-0.31795; -0.30988 ];%für umrechnung der Füllstände in den volt bereich
    k3 = [16395 0 ; 0 15723];% für umrechnung der Pumpleistung in normal Bereich
    k4 = [1.0918 ; 1.0346];%für umrechnung der Pumpleistung in normal Bereich
    D_neu = D(1:100, 1:2);
    C_neu = C(1:100,1:100);
    D_tilde = k3 * Wout * D_neu * k1 * H;
    C_tilde = k3 * Wout * C_neu;
    E = k3 * Wout * D_neu * k1 * k2 + k3*k4;
    disp(size(D));
    D = D_neu;
    C = C_neu;



% Erste Blockreihe
block11 = A_tilde + B_tilde * D_tilde; 
block12 = B_tilde * C_tilde;      

% Zweite Blockreihe
block21 = D * k1 * H;                       
block22 = C;    

% Gesamte Matrix
Matrix = [block11, block12;
                block21, block22];


eigenvals = eigs(Matrix);

magnitudes = abs(eigenvals);
%hier noch einfügen erst ganz normale ridge regression die matrix wird dann
%als start wert benutzt

