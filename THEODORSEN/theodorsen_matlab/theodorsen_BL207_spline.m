clear
clc
close all

%%
Corpo = importXfoilProfile(strcat('BL207.dat'));

x = flipud(Corpo.x);
y = flipud(Corpo.y);

%traslazione del profilo per centrarlo nell'origine
Corpo_traslato = table2array(Corpo);

for i = 1 : length(Corpo.x)
    Corpo_traslato(i, 1) = Corpo_traslato(i, 1)-0.5;
end

%plot del profilo traslato
figure(1)
grid on
h1 = plot(Corpo_traslato(:, 1), Corpo_traslato(:, 2), 'o-b', 'LineWidth', 2);
hold on

%calcolo della linea media col metodo errato, ovvero come media delle
%ordinate di punti non perfettamente allineati sulle ascisse
flag = 0;
y_LM_wrong = [];
for j = 1 : length(Corpo.x)/2    
    y_LM_wrong = [y_LM_wrong; (Corpo_traslato(j, 2) + Corpo_traslato(end-flag, 2))/2];
    flag = flag + 1;
end
x_LM_fwd_wrong = Corpo_traslato(length(Corpo_traslato)/2 + 1 : end, 1);
y_LM_fwd_wrong = flipud(y_LM_wrong);

figure(1)
h2 = plot(x_LM_fwd_wrong, y_LM_fwd_wrong, '-*g');

%% calcolo dell'angolo di Theodorsen
nmax = 4500;
intervalli = 20 : 20 : nmax;

%inizializzazione vettori
alfath_confronto = [];
alfa0_confronto = [];
beta0_confronto = [];

%inizio iterazioni per la valutazione di diverse discretizzazioni della
%linea media, dalla meno fitta alla più fitta. n è il numero di segmenti in
%cui verrà discretizzata
for n = intervalli

    x_spline = linspace(-0.5, 0.5, n);

    y_spline_dorso = cubicspline(flipud(Corpo_traslato(1 : length(Corpo_traslato(:, 1))/2, 1)), flipud(Corpo_traslato(1 : length(Corpo_traslato(:, 1))/2, 2)), x_spline);
    y_spline_ventre = cubicspline(flipud(Corpo_traslato(length(Corpo_traslato(:, 1))/2 + 1 : end, 1)), flipud(Corpo_traslato(length(Corpo_traslato(:, 1))/2 + 1 : end, 2)), x_spline);

    y_LM_fwd_spline = (y_spline_dorso(1:end)+y_spline_ventre(1:end))/2;

    %plot della linea media con 800 punti nella figura con il confronto
    %grafico

        if n == 800
            figure(1)
            h3 = plot(x_spline, y_LM_fwd_spline, '-+r');
        end
        
    I_alfath_vett = [];
    I_alfa0_vett = [];
    I_beta0_vett = [];
    coeff = [];
    a = [];
    b = [];
    b_meno_a = [];
        
        %calcolo degli integrali per alfa_th, alfa_0 e beta_0. Vengono
        %fatte tante iterazioni quanti sono gli intervalli in cui la linea
        %media è stata divisa.
        for i = 1 : length(x_spline) - 1
        
            a = [a; acos(-2*x_spline(i))]; %calcolo degll'estremo di sinistra dell'i-esimo intervalli di integrazione
            b = [b; acos(-2*x_spline(i+1))]; %calcolo degll'estremo di destra dell'i-esimo intervalli di integrazione
            b_meno_a = [b_meno_a; b(i)-a(i)]; %calcolo dell'incremento dell'i-esimo intervallo di integrazione
        
            coeff = [coeff; (y_LM_fwd_spline(i) - y_LM_fwd_spline(i+1))/(x_spline(i) - x_spline(i+1))]; %calcolo della pendenza della linea media (approssimata a tratti lineari) nel tratto i-esimo 

            I = (b(i) - a(i))*(coeff(i)); %calcolo del contributo dell'i-esimo intervallo all'integrale che calcola l'angolo di Theodorsen
            I_alfa0 = ((coeff(i))/pi)*(sin(b(i))-sin(a(i))); %calcolo del contributo dell'i-esimo intervallo all'integrale per il calcolo dell'incidenza con portanza nulla
            I_beta0 = (coeff(i))*((b(i)/2)-(a(i)/2)-0.25*(sin(2*b(i))-sin(2*a(i)))); %calcolo del contributo dell'i-esimo intervallo all'integrale per il calcolo dell'incidenza con momento rispetto all'origine nullo
                
            I_alfath_vett = [I_alfath_vett; I];
            I_alfa0_vett = [I_alfa0_vett; I_alfa0];
            I_beta0_vett = [I_beta0_vett; I_beta0];
         
        end
        
        alfath_rad = (1/pi)*sum(I_alfath_vett); %angolo di Theodorsen, data una linea media discretizzata in n segmenti (radianti)
        alfath_deg = rad2deg(alfath_rad);
        
        alfa0_rad = alfath_rad - sum(I_alfa0_vett); %angolo di incidenza per portanza nulla, data una linea media discretizzata in n segmenti (radianti)        
        alfa0_deg = rad2deg(alfa0_rad); 
        
        beta0_rad = (2/pi)*sum(I_beta0_vett); %angolo di incidenza per momento rispetto all'origine nullo, data una linea media discretizzata in n segmenti (radianti)
        beta0_deg = rad2deg(beta0_rad);
        
        %costruzione vettori per verifica grafica della convergenza del
        %metodo utilizzato
        alfa0_confronto = [alfa0_confronto; alfa0_deg];
        beta0_confronto = [beta0_confronto; beta0_deg];
        alfath_confronto = [alfath_confronto; alfath_deg];

end

figure(1)
legend([h1, h2, h3], 'BL207 GIII', 'Linea media errata', 'Linea media spline', 'Location', 'best')
title('Confronto grafico tra linee medie')
xlabel('x', 'FontWeight', 'bold'); ylabel('y', 'Rotation', 0, 'FontWeight', 'bold')
yline(0); xline(0)
xlim([-0.515, -0.46]); ylim([-0.02, 0.02])
grid on
axis equal

%% scelta del risultato 

toll = 1e-4;
j = 1;
while abs(alfath_confronto(j)-alfath_confronto(j+1)) > toll
    j = j+1;
end
alfa_th = alfath_confronto(j); %Angolo di Theodorsen del profilo

toll = 1e-7;
l = 1;
while abs(alfa0_confronto(l)-alfa0_confronto(l+1)) > toll
    l = l+1;
end
alfa_0 = alfa0_confronto(l); %Angolo di incidenza per portanza del profilo nulla

toll = 1e-6;
m = 1;
while abs(beta0_confronto(m)-beta0_confronto(m+1)) > toll
    m = m+1;
end
beta_0 = beta0_confronto(m); %Angolo di incidenza per momento del profilo rispetto all'origine nullo

fprintf('L''angolo di progetto vale %.5f°\n', alfa_th);
fprintf('L''angolo di zero lift %c%c vale %.4f °\n', char(945), char(8320), alfa_0);
fprintf('L''angolo di momento nullo rispetto all''origine %c%c vale %.4f°\n', char(946), char(8320), beta_0);
%% verifica grafica della convergenza del metodo

%serie di plot per avere una verifica grafica della convergenza del metodo
figure(2)
plot(20:20:nmax, alfath_confronto, '*r')
title('Valori di \alpha_{th} BL207 GIII')
xlabel('n', 'FontWeight', 'bold'); ylabel('\alpha_{th} [°]', 'FontWeight', 'bold')
hold on
plot(intervalli(j), alfa_th, 'y*', 'LineWidth',2)
legend('', '\alpha_{th}',  'Location', 'best')
grid on

figure(3)
plot(20:20:nmax, alfa0_confronto, '*b');
title('Valori di \alpha_{0} BL207 GIII')
xlabel('n', 'FontWeight', 'bold'); ylabel('\alpha_{0} [°]', 'FontWeight', 'bold')
hold on
plot(intervalli(l), alfa_0, 'y*', 'LineWidth',2)
legend('', '\alpha_{0}',  'Location', 'best')
grid on

figure(4)
plot(20:20:nmax, beta0_confronto, '*k');
title('Valori di \beta_{0} BL207 GIII')
xlabel('n', 'FontWeight', 'bold'); ylabel('\beta_{0} [°]', 'FontWeight', 'bold')
hold on
plot(intervalli(m), beta_0, 'y*', 'LineWidth',2)
legend('', '\beta_{0}',  'Location', 'best')
grid on

%% verifica dei risultati (confronto con simulazione XFoil)

CmAC = -(pi/2)*(deg2rad(alfa_0)-deg2rad(beta_0)); %coefficiente di momento riaspetto al centro aerodinamico con incidenza = alfa_0
Cl = 2*pi*(-deg2rad(alfa_0)); %coefficiente di portanza ad incidenza nulla

CmAC_err = 100 - (CmAC*100)/0.0176;
Cl_err = 100 - (Cl*100)/0.1137;

fprintf('Il valore del C_L ad incidenza nulla è %.4f.\n', Cl);
fprintf('Il coefficiente di momento rispetto al centro aerodinamico vale %.5f.\n', CmAC);

%% plot andamento del Cp generato con XFOIL all'incidenza alfa_th
Cp = importXfoilProfile(strcat('cp_theodorsen.dat'));

cp_x = Cp.x;
cp_cp = Cp.y;

figure(5)
plot(cp_x, -cp_cp, '-r', 'LineWidth', 2);
yline(0); xline(0);
xlim([-0.05 1])
ylim([-1 0.6])
pbaspect([3 1 1]);
grid on
% title('Grafico del C_p ad incidenza \alpha_{th}')
xlabel('x', 'FontWeight', 'bold'); ylabel('-C_p', 'FontWeight', 'bold')

saveas(gcf, 'cp_theod.jpg')


