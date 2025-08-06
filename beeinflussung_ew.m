% Fix random seed for comparison purposes
rng(0,'twister');

randVal = rand(1); %sollte 0.8147

%Linearisierung des offenen Systems

%%Arbeitspunkt
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







%% Dreitank Linearisierung

% 0. Parameter
para1.At = 0.0154;
para1.Av = pi*(0.0040)^2;
para1.a12 = 0.3485;
para1.a10 = 0.6677;
para1.g = 9.80665;
para1.a23 = 0.3485;
para1.a30I = 0.6591;
para1.Ts = 0.1;

syms a10 a12 a23 a30I At Av g u1 u2 Ts
x        = sym('x',[3,1]); % Zustandsvektor
y        = sym('y',[3,1]); % Ausgangsvektor
u        = sym('u',[2,1]); % Ausgangsvektor

dx = [(u(1) - a12 * Av * sign(x(1) - x(2)) * sqrt(2*g*abs(x(1) - x(2))) - a10 * Av * sqrt(2*g*x(1))) / At;
    (a12 * Av * sign(x(1) - x(2)) * sqrt(2*g*abs(x(1) - x(2))) - a23 * Av * sign(x(2)-x(3)) * sqrt(2*g*abs(x(2)-x(3)))) / At;
    (u(2) + a23 * Av * sign(x(2)-x(3)) * sqrt(2*g*abs(x(2)-x(3))) - a30I * Av * sqrt(2*g*x(3))) / At];


%% Nichtlineares Modell (als Funktion)
dreitank_model = @(x, u, para) [
    (u(1) - para1.a12 * para1.Av * sign(x(1) - x(2)) * sqrt(2*para1.g*abs(x(1) - x(2))) ...
     - para1.a10 * para1.Av * sqrt(2*para1.g*x(1))) / para1.At;

    (para1.a12 * para1.Av * sign(x(1) - x(2)) * sqrt(2*para1.g*abs(x(1) - x(2))) ...
     - para1.a23 * para1.Av * sign(x(2)-x(3)) * sqrt(2*para1.g*abs(x(2)-x(3)))) / para1.At;

    (u(2) + para1.a23 * para1.Av * sign(x(2)-x(3)) * sqrt(2*para1.g*abs(x(2)-x(3))) ...
     - para1.a30I * para1.Av * sqrt(2*para1.g*x(3))) / para1.At
];

%% Diskrete Übergangsfunktion (Euler/ODE45)
simulate_step = @(x, u, T, para1) ...
    ode45_step(@(t, x_) dreitank_model(x_, u, para1), x, T);

function x_next = ode45_step(f, x0, T)
    [~, x_out] = ode45(f, [0 T], x0);
    x_next = x_out(end,:)';
end




%% Neuer Arbeitspunkt: numerische Bestimmung von xR mit f_d(xR, uR) = xR

% Funktion f_d wie vorher (diskretisiertes nichtlineares Modell)
f_d = @(x, u) simulate_step(x, u, T, para1);  % simulate_step wie vorher definiert

% Wir suchen xR: x = f_d(x, uR) → also: f_fix(x) = f_d(x, uR) - x = 0
f_fix = @(x) f_d(x, uR) - x;

% Startwert z. B. dein alter xR
xR_guess = xR;

% fsolve zur Bestimmung des Fixpunkts
options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-10);
[xR_new, fval, exitflag] = fsolve(f_fix, xR_guess, options);

 % if exitflag <= 0
 %     warning('fsolve hat keinen gültigen Arbeitspunkt gefunden.');
 % else
 %     disp('Neuer Arbeitspunkt xR:');
 %     disp(xR_new);
 % end

% Setze neuen Arbeitspunkt
xR = xR_new;
vxr = [xR(1);xR(3);xR(1);xR(3)];

 %% Numerische Linearisierung um Arbeitspunkt
 delta = 1e-6;
 n = length(xR); m = length(uR);
 A_d = zeros(n);
 B_d = zeros(n, m);
 
 % A-Matrix
 for i = 1:n
     dx = zeros(n,1); dx(i) = delta;
     f1 = f_d(xR + dx, uR);
     f2 = f_d(xR - dx, uR);
     A_d(:,i) = (f1 - f2) / (2*delta);
 end
 
 % B-Matrix
 for i = 1:m
     du = zeros(m,1); du(i) = delta;

     f1 = f_d(xR, uR + du);
     f2 = f_d(xR, uR - du);
     B_d(:,i) = (f1 - f2) / (2*delta);
 end









%laden der Daten und splitten der Daten
data = readmatrix('daten_simulation_v13.dat');%testen mit v4 und v5, v6 erzielt beim anlernen keine besseren ergebnisse und ist deutlich rechenaufwändiger, da dort 80.000 datenpunkte drinnenliegen

x1 = data(:,3);
x3 = data(:,4);
u1 = data(:,1);
u3 = data(:,2);

vx1 = (x1-0.31795)/-0.3314;
vx3 = (x3-0.30988)/-0.33781;

vu1 = 16395*u1-1.0918;
vu3 = 15723*u3-1.0346;


% Clipping (Begrenzen auf [-1, 1])
vx1 = min(max(vx1, -1), 1);
vx3 = min(max(vx3, -1), 1);
vu1 = min(max(vu1, -1), 1);
vu3 = min(max(vu3, -1), 1);

x = [vx1 , vx3];
u = [vu1 , vu3];

t = 1;%delay time 

X_t = x(1:end-t,:);%aktuell gewähltes Delay 1
X_t_dt = x(1+t:end, :);

%X_dot_t = (X_t - X_t_dt)/t;
U_t = u(1:end-t,:);

eingang = [X_t ,X_t_dt]';
target = U_t';





%erstellen des ESN und festlegen der Parameter
esn = ESN();
esn.inputs = 4;
esn.input_scaling = 0.85;
esn.neuronen = 100;
esn.state = zeros(esn.neuronen,1);
esn.conectivity = 0.3;
esn.reg = 1e-4;
esn.leak_rate = 0.75;
esn.spectral_radius = 0.95;
esn.washout = 100;
esn.use_washout = true;
esn.use_bias = true;
esn.bias = 0.8;
esn.use_bias_out = true;
esn.bias_out = 1.0;%bei tests hier keinen einfluss auf den mse gehabt => es würde reichen nur einen input bias zu haben
esn.initialize();


eingang1 = vxr .* ones(4,50000);
um = [vxr;esn.bias*esn.input_scaling];
states = esn.run(eingang1);
statem = states(:,49900);
xm = statem;

[C,D] = linearizeESN(esn.Wres,esn.Win,xm,um,esn.leak_rate);
% hier regelkreis ohne referenz oder bias genutzt, weil diese eh nicht
% beeinflussbar
H = [1 0 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
k1 = [1/-0.3314 0 0 0; 0 1/-0.33781 0 0; 0 0 1/-0.3314 0; 0 0 0 1/-0.33781];
k2 = [-0.31795; -0.30988;-0.31795; -0.30988];
k3 = [1/16395 0 ; 0 1/15723];% für umrechnung der Pumpleistung in normal Bereich
k4 = [1.0918 ; 1.0346];%für umrechnung der Pumpleistung in normal Bereich
D_neu = D(1:100, 1:4);
C_neu = C(1:100,1:100);



% Zweite Blockreihe
block21 = D_neu * k1 * H;                       
block22 = C_neu;    

states_run = esn.run(eingang);
esn.train(target,states_run);
states_run = [states_run;esn.bias*ones(1, size(states_run,2))];
target = target(:,esn.washout+1:end);

% opts = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',1e4);

opts = optimoptions('lsqnonlin', ...
    'Algorithm','interior-point',...
    'Display', 'iter', ...
    'FunctionTolerance', 1e-8, ...
    'StepTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-8, ...
    'MaxFunctionEvaluations', 1e5);


A_tilde = [A_d(1,1) A_d(2,1) A_d(3,1)  0 0 ; A_d(1,2) A_d(2,2) A_d(3,2)  0 0 ; A_d(1,3) A_d(2,3) A_d(3,3)  0 0 ; 0 0 0 1 0 ; 0 0 0 0 1];
B_tilde = [B_d;0 0;0 0];


para.H = H;
para.D_neu = D_neu;
para.C_neu = C_neu;
para.k1 = k1;
para.k2 = k2;
para.k3 = k3;
para.k4 = k4;
para.block21 = block21;
para.block22 = block22;
para.n_out = size(target,1);
para.n_states = (size(states_run,1));
para.r = 0.8;
para.A_tilde = A_tilde;
para.B_tilde = B_tilde;
washout = esn.washout;
ridge = esn.reg;
Wout = 0*esn.Wout;
Wout_vec = Wout(:);





function err = ridgeregression(Wout_vec,states,ziel,ridge) % zu optimierende funktion für lsqnonlin
    n_out = size(ziel,1);
    n_states = (size(states,1));
    
    Wout = reshape(Wout_vec, n_out,n_states);

    fehler = ziel - Wout * states;
    fitting = sqrt(ridge) * Wout;
    err = [fehler(:);fitting(:)];
end

 function [c,ceq] = nebenbed(Wout_vec,params)
     n_out = params.n_out;
     n_states = params.n_states;
     r = params.r;
     Wout = reshape(Wout_vec,n_out,n_states);
     Wout = Wout(1:2,1:100);
 
     %Matrizen aktualisieren
     D_tilde = params.k3 * Wout * params.D_neu * params.H;
     C_tilde = params.k3 * Wout * params.C_neu;
     %Blöcke erstellen/aktualisieren
     block11 = params.A_tilde + params.B_tilde*D_tilde;
     block12 = params.B_tilde * C_tilde;
     block21 = params.block21;
     block22 = params.block22;
     %aufstellen matrix
     M = [block11,block12;
         block21,block22];
     eigenvalues = eig(M);
     magnitudes = abs(eigenvalues);
 
     c1 = magnitudes - r;
     c = c1(4:end-2);
     ceq = [];
 end
 
 %aufruf von lsqnonlin
 Wopt_vec = lsqnonlin( ...
     @(W) ridgeregression(W, states_run, target, ridge), ...
     Wout_vec, ...
     [], [], [],[],[],[],...
     @(W) nebenbed(W, para),...
      opts);
     
    
 Wopt = reshape(Wopt_vec,para.n_out,para.n_states);
esn.use_washout = false;
esn.Wout = Wopt;
esn.reset_state();
out_test = esn.compute_output(eingang);
out_test = out_test(:,washout+1:end);

%berechnen der Fehler zum vergleichen
mse = esn.MSE(target,out_test);

fehler_durchschnitt = mean(mse,'all');

fehler_summe = sum(mse,'all');


%Matrizen aktualisieren
Wopt = Wopt(1:2,1:100);
     D_tilde = para.k3 * Wopt * para.D_neu * para.H;
     C_tilde = para.k3 * Wopt * para.C_neu;
     %Blöcke erstellen/aktualisieren
     block11 = para.A_tilde + para.B_tilde*D_tilde;
     block12 = para.B_tilde * C_tilde;
     block21 = para.block21;
     block22 = para.block22;
     %aufstellen matrix
     M = [block11,block12;
         block21,block22];
     eigenvalues = eig(M);
     magnitudes = abs(eigenvalues);
    

    %  X = [];
    % D_neu = D(:,1:4);
    % for k = 1:esn.neuronen -1
    %     temp = C_neu^k * D_neu;
    %     X = [X, temp];
    % 
    % 
    % end
