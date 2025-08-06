function [A_lin, B_lin] = linearizeESN(Wres, Win, xm, um, alpha)%hier ist alpha die Leakrate, Win die Eingangsmatris, Wre die Reservoirmatrix, xm die internen Zustände des ESN in ruhe, um die System zustände in ruhe
    % Anzahl der Neuronen
    N = size(Wres, 1);
    K = length(um);  % Anzahl der Eingänge

    % Initialisieren
    A_lin = zeros(N, N);
    B_lin = zeros(N, K);  % Jetzt Matrix, nicht Vektor
    
    %Berechnung folgt der in der Bachelorarbeit angegebenen Quelle
    % Berechne z_m = Win * um + Wres * xm
    z_m = Win * um + Wres * xm;  % Größe N x 1

    % Berechne sech^2(z) = 1 - tanh(z)^2
    sech2 = 1 - tanh(z_m).^2;  % Elementweise

    % Baue A_lin (N x N)
    for i = 1:N
        A_lin(i, :) = Wres(i, :) * sech2(i);
    end
    A_lin = (1 - alpha) * eye(N) + alpha * A_lin;

    % Baue B_lin (N x K)
    for i = 1:N
        for k = 1:K
            B_lin(i, k) = Win(i, k) * sech2(i);
        end
    end
    B_lin = alpha * B_lin;
end
