# Bachelorarbeit-Reservoir-Computing-basierte-Regelung-linearer-und-nichtlinearer-Systeme
Dieses Repository enthält alle für meine Bachelorarbeit notwendigen Codeabschnitte, darunter eine Matlab-Klasse zur Erstellung von ESNs.


- ESN.m file ist die Grundfile die die Methoden zur Erstellung eines ESN beinhaltet (done)
- Stabilitaets_test.m enthält den Code der für die Stabilitätstestung des geschlossenen Regelkreises verwendet werden kann hier wird in Abhängigkeiten ob Biase verwendet werden die Dynamikmatrix des geschlossenen Regelkreises berechnet
- stabilitaetstest_ohneref.m enthält den Code der für die Stabilitätsrechnung tatsächlich genutzt wurde hier werden weder Referenz noch Bias mit in die Dynamikmatrix mit einbezogen
- linearizeESN.m beinhaltet den Code zur Linearisierung des ESN dies wird für den Stabilitätstest und die Beeinflussung der Eigenwerte benötigt
- nutzung_ESN.m ist der Code der für das Erstellen der ESN genutzt wurde (done)
- Trainingsdaten erstellen enthält den Code der genutzt wurde zur Erstellung der hochwertigen Trainingsdaten (done)
- die hochwertigen Trainigsdaten befinden sich in daten_simulation_v13.dat, die fehlerhaften Trainingsdaten in daten_simulation_v5 (done)
- eigenwerte_ohnenb.m ist der Code zur Berechnung der Ausgangsmatrix während man die Eigenwerte des geschlossenen Regelkreises beeinflusst

