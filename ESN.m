classdef ESN < handle
    properties
        neuronen %anzahl der neuronen default 100
        inputs %anzahl der inputs default 4
        outputs %anzahl der outputs default 2
        conectivity %verbindungen in dem reservoir default0.2
        input_scaling %default 1.0
        leak_rate %default 0.3 gibt an wie schnell das Modell alte zustände vergisst bei 1 linear=> sehr schnell
        spectral_radius %größter EW von self.Wres default 0.9
        bias %default 1.0 vermeidung von overfitting, für input
        bias_out %default 1.0 für output
        Win %Matrix von inputs ins reservoir
        Wout %output Mmatrix
        Wres % reservoir matrix
        reg %wert für ridge regression default 1e-6
        state %aktueller zustand der neuronen default startwert ist 0
        washout;%anzahl an ignorierten Datenpunkten default 50
        use_washout;%ob washput genutzt wird defalut true
        use_bias;%ob bias genutzt werden sollen bei input
        use_bias_out;%ob bias auch in state zu output genutzt wird
    end
    methods(Access = public)
        
        function self = default(self)%function die alle werte im default initialisiert
            self.neuronen = 100;
            self.inputs = 4;
            self.outputs = 2;
            self.conectivity = 0.2;
            self.input_scaling = 1.0;
            self.leak_rate = 0.3;
            self.spectral_radius = 0.9;
            self.bias = 1.0;
            self.use_bias = true;
            self.use_bias_out = false;
            self.Win = [];
            self.Wres = [];
            self.Wout = [];
            self.reg = 1e-6;
            self.state = zeros(self.neuronen,1);
            self.washout = 10;
            self.use_washout = true;
            
        end
        
        function self = ESN()%aufrufen um alles mit default zu initialisieren
            self.default();
            %alles auf den default setzen


        end

        function initialize(self)%muss aufgerufen werden um die Matrizen Wres und Win zu berechnen
            if self.use_bias
                self.Win = (rand(self.neuronen,self.inputs+1)-0.5)*2*self.input_scaling;%initialisierung von Win
            else
                self.Win = (rand(self.neuronen,self.inputs)-0.5)*2*self.input_scaling;%initialisierung von Win
            end
            self.Wres = sprand(self.neuronen,self.neuronen,self.conectivity);
            self.Wres(self.Wres~= 0) = self.Wres(self.Wres ~=0)-0.5;
            self.Wres = full(self.Wres);
            opts.tol = 1e-3;
            eigenvals = eigs(self.Wres,1,'LM',opts);
            self.Wres = self.Wres*(self.spectral_radius/abs(eigenvals));



        end
        function test_stabil(self)%betrachtung des spectralradius der effektiven matrix
            h = eigs((1-self.leak_rate)*eye(self.neuronen)+self.leak_rate*self.Wres);
            if h < 1
                disp("Der Spectralradius der effektiven Matrix ist kleiner als 1, das ESN ist stabil");
            else
                disp("Der Spectralradius der effektiven Matrix ist größer als 1, das ESN ist nicht stabil");
            end


        end
        function test_steuerbar(self)%hier wird die Steuerbarkeit anhand des linearen ESN getestet wie hier implementiert wird um x = 0; u= 0 getestet
            
            % Leaky Rate
            alpha = self.leak_rate;  
        
            % Effektive Dynamikmatrix
            A_eff = (1 - alpha) * eye(self.neuronen) + alpha * self.Wres;
        
            % Initialisierung der Steuerbarkeitsmatrix
            C = zeros(self.neuronen, self.neuronen * size(self.Win, 2));
            C(:, 1:size(self.Win, 2)) = self.Win;
        
            % Konstruktion der Steuerbarkeitsmatrix mit A_eff
            for k = 1:self.neuronen - 1
                start_idx = k * size(self.Win, 2) + 1;
                end_idx = (k + 1) * size(self.Win, 2);
                C(:, start_idx:end_idx) = A_eff * C(:, (k - 1) * size(self.Win, 2) + 1 : k * size(self.Win, 2));
            end
        
            % Berechnung des Rangs
            rank_C = rank(C);
            
            % Ausgabe
            if rank_C == self.neuronen
                disp("Das System ist vollständig steuerbar in der Nähe von x = 0");
            else
                disp("Das System ist NICHT vollständig steuerbar in der nähe von x = 0");
            end
        end

        function states = run(self,eingang)%berechnet aktuelle Zustände im Reservoir
            T = length(eingang);
            x_n = self.state;%initialer zustand
            X = zeros(self.neuronen, T);
            
            for t = 1:T
                x_t = eingang(:,t);
                if self.use_bias
                    x_bias = [x_t ; self.bias];
                else
                    x_bias = x_t;
                end
                x_n = (1-self.leak_rate)*x_n + self.leak_rate* tanh(self.Wres*x_n + self.Win *x_bias);%neuronen zustände berechnen
                X(:,t) = x_n;
            end
            self.state = x_n;
            if self.use_washout
                states = X(:,self.washout+1:end);
            else
                states = X;
            end

        end
        function self = train(self,ziel,states)%trainiert die readout matrix mit ridge regression(analytisch)
            if self.use_washout
                target = ziel(:,self.washout+1:end);
            else
                target = ziel;
            end
            if self.use_bias_out
                X_aug = [states;self.bias_out*ones(1,size(states,2))];
            else
                X_aug = states;
            end
            self.Wout = target * X_aug'/(X_aug * X_aug' + self.reg * eye(size(X_aug,1)));

        end
        function output = compute_output(self,eingang)%erstellt den ausgang des ESN zu einem Eingang/nach dem training des ESN
            states = self.run(eingang);
            
            if self.use_bias_out
                X_aug = [states;self.bias*ones(1,size(states,2))];
            else
                X_aug = states;
            end
            output = self.Wout*X_aug;
        end
        function reset_state(self)%resetet die internen Zustände der Neuronen auf 0
            self.state = zeros(self.neuronen,1);
        end
        function fehler = MSE(self,out,output)%berechnet den MSE
                if self.use_washout
                    target = out(:,self.washout+1:end);
                else
                    target = out;
                end
                    fehler = mean((target-output).^2);
        end
        function self = complete_train(self,eingang,target)%etwas vereinfachtere Trainingsvariante/direkt vollständiges Training man erhält allerdings keine internen States zurück
            states = self.run(eingang);
            self.train(target,states);
        end
    end
    

end