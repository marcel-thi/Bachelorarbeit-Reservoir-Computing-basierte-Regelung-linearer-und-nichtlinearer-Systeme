% Fix random seed for comparison purposes
rng(0,'twister');

randVal = rand(1); %sollte 0.8147

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

t = 1;%delay time im Realsystem auf 15 gestzt

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
esn.bias_out = 0.8;%bei tests hier keinen einfluss auf den mse gehabt => es würde reichen nur einen input bias zu haben
esn.initialize();

states = esn.run(eingang);

esn.train(target,states);
%esn.complete_train(eingang,target);
%esn.reset_state();%setzt state der neuronen auf 0 zurück
load('last_state_ESN.mat');
esn.state = last_state_ESN; %ein state der neuronen aus der simulation mit washout kein Unterschied ohne großer unterschied

%ausprobieren an test werten
esn.use_washout = true;
esn.washout = 100;
out_test = esn.compute_output(eingang);


%berechnen der Fehler zum vergleichen
mse = esn.MSE(target,out_test);

fehler_durchschnitt = mean(mse,'all');

fehler_summe = sum(mse,'all');

%speichern aller Parameter die in simulink gebraucht werden
p.Wout = esn.Wout;
p.Win = esn.Win;
p.Wres = esn.Wres;
p.leak_rate = esn.leak_rate;
%p.state = esn.state;
p.neuronen = esn.neuronen;

p.bias = esn.bias;
p.bias_out = esn.bias_out;
p.state = last_state_ESN;

params = Simulink.Parameter;
params.Value = p;


