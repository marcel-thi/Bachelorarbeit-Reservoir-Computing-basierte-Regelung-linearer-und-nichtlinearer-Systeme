
%%Arbeitspunkt
xr1 =  0.25;
xr2 =  0.2;
xr3 =  0.15;
xR = [xr1;xr2;xr3];
ur1 =  0.0001;
ur3 =  0.0;
uR = [ur1;ur3];

vxr1 = (xr1-0.31795)/-0.3314;
vxr3 = (xr3-0.30988)/-0.33781;
vxr = [vxr1;vxr3;vxr1;vxr3];






%% Dreitank Linearisierung

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
y        = sym('y',[2,1]); % Ausgangsvektor
u        = sym('y',[2,1]); % Ausgangsvektor

dx = [(u(1) - a12 * Av * sign(x(1) - x(2)) * sqrt(2*g*abs(x(1) - x(2))) - a10 * Av * sqrt(2*g*x(1))) / At;
    (a12 * Av * sign(x(1) - x(2)) * sqrt(2*g*abs(x(1) - x(2))) - a23 * Av * sign(x(2)-x(3)) * sqrt(2*g*abs(x(2)-x(3)))) / At;
    (u(2) + a23 * Av * sign(x(2)-x(3)) * sqrt(2*g*abs(x(2)-x(3))) - a30I * Av * sqrt(2*g*x(3))) / At];

x_next = x + Ts * dx;

% Systemmatrizen
A = jacobian(x_next,x);
b = jacobian(x_next,u);

A = subs(A,x,xR); 
b = subs(b,x,xR);

A = subs(A,u,uR); 
b = subs(b,u,uR);

para.u1R = ur1; % Definition des Arbeitspunkts, muss physikalisch sinnvoll sein 
para.u2R = ur3;
para_save = para; 
syms_params = {At, Av, a12, a10, a23, a30I, g,Ts};
para_values = {para.At, para.Av, para.a12, para.a10, para.a23, para.a30I, para.g,para.Ts};

AN = double(subs(A, syms_params, para_values));
bN = double(subs(b, syms_params, para_values));

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
Wout = esn.Wout;
if esn.use_bias && esn.use_bias_out == false
    disp("mit input bias")
    H = [1 0 0 0 0 0 ; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0 ; 0 0 0 0 0 1]; 
    A_tilde = [AN(1,1) AN(2,1) AN(3,1)  0 0 0; AN(1,2) AN(2,2) AN(3,2)  0 0 0; AN(1,3) AN(2,3) AN(3,3)  0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
    B_tilde = [bN ; 0 0; 0 0; 0 0];
    k1 = [1/-0.3314 0 0 0 0; 0 1/-0.33781 0 0 0; 0 0 1/-0.3314 0 0; 0 0 0 1/-0.33781 0; 0 0 0 0 1];% für umrechnung der Füllstände in volt bereich
    k2 = [-0.31795; -0.30988 ; -0.31795; -0.30988;0];%für umrechnung der Füllstände in den volt bereich
    k3 = [16395 0 ; 0 15723];% für umrechnung der Pumpleistung in normal Bereich
    k4 = [1.0918 ; 1.0346];%für umrechnung der Pumpleistung in normal Bereich
    D_tilde = k3 * Wout * D * k1 * H;
    C_tilde = k3 * Wout * C;
    E = k3 * Wout * D * k1 * k2 + k3*k4;
elseif esn.use_bias && esn.use_bias_out
    H = [1 0 0 0 0 0 ; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0 ; 0 0 0 0 0 1]; 
    A_tilde = [AN(1,1) AN(2,1) AN(3,1)  0 0 0; AN(1,2) AN(2,2) AN(3,2)  0 0 0; AN(1,3) AN(2,3) AN(3,3)  0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0 ; 0 0 0 0 0 1];
    B_tilde = [bn;0 0 ;0 0 ; 0 0];
    k1 = [1/-0.3314 0 0 0 0; 0 1/-0.33781 0 0 0; 0 0 1/-0.3314 0 0; 0 0 0 1/-0.33781 0; 0 0 0 0 1];% für umrechnung der Füllstände in volt bereich
    k2 = [-0.31795; -0.30988 ; -0.31795; -0.30988;0];%für umrechnung der Füllstände in den volt bereich
    k3 = [16395 0 ; 0 15723];% für umrechnung der Pumpleistung in normal Bereich
    k4 = [1.0918 ; 1.0346];%für umrechnung der Pumpleistung in normal Bereich
    D_neu = [D;zeros(1,size(D,2))];
    C_neu = [C; zeros(1,size(C,2))];
    C_neu = [C_neu, zeros(1,size(C_neu,1))'];
    C_neu(size(C_neu,1),size(C_neu,2)) = 1;
    D_tilde = k3 * Wout * D_neu * k1 * H;
    C_tilde = k3 * Wout * C_neu;
    E = k3 * Wout * D_neu * k1 * k2 + k3*k4;
    D = D_neu;
    C = C_neu;
else
    H = [1 0 0 0 0 ;0 0 1 0 0 ; 0 0 0 1 0; 0 0 0 0 1];
    disp("ohne inpiut bias")
    A_tilde = [AN(1,1) AN(2,1) AN(3,1)  0 0 ; AN(1,2) AN(2,2) AN(3,2)  0 0 ; AN(1,3) AN(2,3) AN(3,3)  0 0 ; 0 0 0 1 0; 0 0 0 0 1];
    B_tilde = [bN ; 0 0; 0 0];
    k1 = [1/-0.3314 0 0 0 ; 0 1/-0.33781 0 0 ; 0 0 1/-0.3314 0 ; 0 0 0 1/-0.33781 ];% für umrechnung der Füllstände in volt bereich
    k2 = [-0.31795; -0.30988 ; -0.31795; -0.30988];%für umrechnung der Füllstände in den volt bereich
    k3 = [16395 0 ; 0 15723];% für umrechnung der Pumpleistung in normal Bereich
    k4 = [1.0918 ; 1.0346];%für umrechnung der Pumpleistung in normal Bereich
    D_tilde = k3 * Wout * D * k1 * H;
    C_tilde = k3 * Wout * C;
    E = k3 * Wout * D * k1 * k2 + k3*k4;
end




% Erste Blockreihe
block11 = A_tilde + B_tilde * D_tilde;  % Größe z.B. 5x5
block12 = B_tilde * C_tilde;            % Größe z.B. 5x2

% Zweite Blockreihe
block21 = D * k1 * H;                       % z.B. 2x5
block22 = C;                            % z.B. 2x2

% Gesamte Matrix
Matrix = [block11, block12;
                block21, block22];


eigenvals = eigs(Matrix);

magnitudes = abs(eigenvals);


