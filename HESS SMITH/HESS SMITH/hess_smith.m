%% traccia Hess Smith (2024)

clc
close all
clear 

% addpath mat_functions

%% 
% Questo script implementa il metodo di Hess Smith come fornito dalla
% consegna, successivamente alla risoluzione del sistema lineare sono state
% aggiunte le righe necessarie a determinare il campo di velocità attorno
% ad un profilo NACA 0012. 
% Una volta determinato il campo di velocità si è determinato il campo di
% pressione tramite il teorema di Bernoulli.
% Si fa anche un confronto tra i dati ottenuti con questo script e i dati
% provenienti da un implementazione in xfoil dello stesso profilo. 
% I calcoli vengono effettuati per un angolo di incidenza di 2° e per una
% velocità della corrente incidente di 1 m/s.
% La pannellizzazione presenta 101 pannelli, il profilo viene costruito con
% corda unitaria e leading edge centrato nell'origine. 


%% Input

U_inf = 1;  % Velocità all'infinito [m/s]
alpha = 2 ;   % Angolo di incidenza [°]
U_inf_x = U_inf * cos(deg2rad(alpha));
U_inf_y = U_inf * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y];
U_inf_normal = [-U_inf(2); U_inf(1)];
U_inf_normal = U_inf_normal ./ norm(U_inf_normal);

TestCase = 0;

CodiceProfilo = '0012';
Chord = 1;
NPannelli = 101;

LE_X_Position = 0;
LE_Y_Position = 0;

%% Creazione profilo

% numero profilo:
% [x,y]=createProfile(CodiceProfilo,NPannelli,Chord);

Corpo = importXfoilProfile(strcat('NACA_', CodiceProfilo, '.dat'));
% Prima flippa i vettori
x = flipud(Corpo.x);
y = flipud(Corpo.y);
Corpo.x = x.*Chord;
Corpo.y = y.*Chord;

figure (1);
plot(x, y, '-o r', 'LineWidth',0.8)
grid on
title('Rappresentazione Profilo');
xlabel('Corda', 'FontWeight','bold')
axis equal
saveas(gcf, 'rappresentazione_profilo.jpg')
%% Creazione di una struttura di pannelli

[Centro, Normale, Tangente, Estremo_1, Estremo_2, alpha_1, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);
        
%% Inizializzazione matrici e vettori

% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori

NCols = sum(NPannelli) + 1;
NRows = NCols;
matriceA = zeros(NRows, NCols);
TermineNoto = zeros(NRows, 1);

%% Creazione della matrice quadrata As


for i = 1:NPannelli
    index_i = i; % riga

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;  % Colonna

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            matriceA(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);

            matriceA(index_i, sum(NPannelli)+1) = matriceA(index_i, sum(NPannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);


        end

end


%% Creazione delle componenti dei vettori a_v, c_s e c_v


Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';


b = 0;
for j = 1:NPannelli(1)

    index_j = j;

    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    b = b + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    b = b + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);


    matriceA(sum(NPannelli) + 1, index_j) = a;

end

matriceA(sum(NPannelli) + 1, sum(NPannelli) + 1) = b;



%% Creazione del termine noto

for j = 1:NPannelli

    Normale_qui = Normale(j, :)';

    index = j;

    TermineNoto(index) = - dot(U_inf, Normale_qui);
end

Tangente_1 = Tangente(1, :)';
Tangente_end = Tangente(end, :)';
TermineNoto(sum(NPannelli) + 1) = - dot(U_inf, (Tangente_1 + Tangente_end));

%% Risoluzione sistema lineare
Soluzione = linsolve(matriceA,TermineNoto);


%% Calcolo del campo di velocità attorno al profilo 
[Centro, Normale, Tangente, Estremo_1, Estremo_2, alpha, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);

velU_sorgenti = zeros(NPannelli,2);
velU_vortici = zeros(NPannelli,2);
U_TOT = zeros(NPannelli,2);

for i = 1:NPannelli
    index_i = i; % riga

    Centro_qui = Centro(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;  % Colonna

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            velU_sorgenti(index_j,:) = Soluzione(j).*ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
            velU_vortici(index_j,:) = Soluzione(end).*ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);

        end
            

    U_TOT(i,:) = U_inf' + sum(velU_vortici) + sum(velU_sorgenti);
end

% campo di velocità attorno al profilo, la velocità è stata calcolata nel
% centro di ogni pannello. 


%% Calcolo del campo di pressione 

%dati necessari per usare il teorema di bernoulli, la corrente incidente è
%considerata in condizioni standard 

P_inf = 101325; %[Pa]
rho_inf = 1.225; %[kg/m^3]
P = zeros(NPannelli,1); 
cp = zeros(NPannelli,1); 

for i = 1 : NPannelli
    index_i = i;
    P(i) = P_inf + (0.5*rho_inf) - 0.5*rho_inf* (norm(U_TOT(i,:))^2);
    cp(i) = (P(i) - P_inf)./(0.5*rho_inf);


end

% Calcolo dei coefficienti di pressione, suddiviso tra dorso e ventre: 

cp_dorso = -cp(51:end,1);
x_dorso = x(52:end,1);

cp_ventre = -cp(1:51,1);
x_ventre = x(1:51);


%% Confronto con i risultati ottenuti da XFOIL per il NACA0012

cp_xfoil = load("CPWR_NACA_0012_2.txt"); 
% questo file contiene i dati del grafico del cp in corda nelle stesse
% condizioni analizzate in questo script, ottenuto tramite xfoil 


figure (2)
plot( x_dorso, cp_dorso, '*-', 'LineWidth',0.8)
grid on
hold on
plot( x_ventre, cp_ventre, '*-','LineWidth',0.8)
plot(cp_xfoil(1:51,1),-cp_xfoil(1:51,2),'LineWidth',0.8)
plot(cp_xfoil(52:end,1),-cp_xfoil(52:end,2),'LineWidth',0.8)
legend ('C_P dorso Matlab','C_P ventre Matlab','C_P dorso Xfoil','C_P ventre Xfoil');
% title ('Confronto Cp Xfoil vs Matlab');
xlabel('x', 'FontWeight','bold'); ylabel('-C_P', 'FontWeight','bold');

saveas(gcf, 'confronto_cp.jpg')


%% Calcolo del CL con il teorema di Bernoulli 

cl_bernoulli = 0;
U_inf_normal = [-U_inf(2); U_inf(1)];
cpp = -cp;
for i = 1 : NPannelli

    cl_bernoulli = cl_bernoulli  + cpp(i)*lunghezza(i)*(dot(Normale(i, :), U_inf_normal)); 

end



%% Calcolo del CL con il teorema di Kutta Joukowsky

Gamma = 0;

for i = 1 : NPannelli
    Gamma = Gamma + lunghezza(i)*(Soluzione(end));
end

cl_KJ = 2*Gamma; 


%% confronto con xfoil del cl 
data_cl = load('polar_cl.txt');
cl_xfoil = data_cl(5,2);

fprintf('Il coefficiente di portanza ad alpha 2° ottenuto con il teorema di Bernoulli è %f\n' , cl_bernoulli);

fprintf('Il coefficiente di portanza ad alpha 2° ottenuto con il teorema di Kutta Joukowsky è %f\n' , cl_KJ);

fprintf('Il coefficiente di portanza ad alpha 2° ottenuto da xfoil è %f\n' , cl_xfoil);

%%  Calcolo del coefficiente di momento rispetto al bordo di attacco e rispetto al centro aerodinamico 
cm_LE = 0;
CentroCM = [Centro, zeros(NPannelli, 1)];
NormaleCM = [Normale, zeros(NPannelli, 1)];
z = [0,0,1]';

for i = 1 : NPannelli
 
     cm_LE = cm_LE + cpp(i)*lunghezza(i)*cross(CentroCM(i,:),NormaleCM(i,:));
end

cm_LE = -cm_LE(3);

% il centro aerodinamico viene considerato ad 1/4 della corda
cm_AC = cm_LE + cl_KJ *0.25; 

%% confronto con xfoil 

data_cm = load('polar_cm.txt');
cm_xfoil = data_cm(5,2);
fprintf('Il coefficiente di momento ad alpha 2° ottenuto da matlab è %f\n' , cm_AC);
fprintf('Il coefficiente di portanza ad alpha 2° ottenuto da xfoil è %f\n' , cm_xfoil);
