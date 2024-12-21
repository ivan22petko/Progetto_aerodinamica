%% traccia Hess Smith (2024)

clc
close all
clear 

% addpath mat_functions

% in questo script il metodo di hs è implementato per una serie di valori
% di incidenza alpha variabili tra 0 e 5 gradi, in modo tale da potere
% costruire la curva CL/alpha. 

%% Input
alpha = (0:0.5:5) ;   % Angolo di incidenza [°]
CM=[];
CL=[];

%% 
for e=1:length(alpha)


    U_inf = 1;  % Velocità all'infinito [m/s]
    U_inf_x = U_inf * cos(deg2rad(alpha(e)));
    U_inf_y = U_inf * sin(deg2rad(alpha(e)));
    
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
    
    %% Calcolo del cp e della velocità sui pannelli
    [Centro, Normale, Tangente, Estremo_1, Estremo_2, ~, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);
    
    velU_sorgenti = zeros(NPannelli,2);
    velU_vortici = zeros(NPannelli,2);
    U_TOT = zeros(NPannelli,2);
    
    for i = 1:NPannelli
        index_i = i; % riga
        
        % posizione in cui mi voglio mettere ovvero in cui calcolo il campo di
        % velocità 
    
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
    
    P_inf = 101325; 
    rho_inf = 1.225; 
    P = zeros(NPannelli,1); 
    cp = zeros(NPannelli,1); 
    
    for i = 1 : NPannelli
        index_i = i;
        P(i) = P_inf + (0.5*rho_inf) - 0.5*rho_inf* (norm(U_TOT(i,:))^2);
        cp(i) = (P(i) - P_inf)./(0.5*rho_inf);
    
    
    end
       
    
    %% Calcolo del CL 
    cl1 = 0;
    U_inf_normal = [-U_inf(2); U_inf(1)];
    cpp = -cp;
    for i = 1 : NPannelli
    
        cl1 = cl1 + cpp(i)*lunghezza(i)*(dot(Normale(i, :), U_inf_normal)); %modulo giusto. Segno NO
        CL(e)=cl1;

    end
      
    
    %%  Calcolo del CM
    cm_LE = 0;
    CentroCM = [Centro, zeros(NPannelli, 1)];
    NormaleCM = [Normale, zeros(NPannelli, 1)];
    z = [0,0,1]';
    
    for i = 1 : NPannelli
     
         cm_LE = cm_LE + cpp(i)*lunghezza(i)*cross(CentroCM(i,:),NormaleCM(i,:));
    end
    
    cm_LE = -cm_LE(3);
    
    cm_AC = cm_LE + cl1*0.25;

    CM(e)=cm_AC;

end


%% Confronto grafico del CL 

data_cl = load('polar_cl.txt'); %dati presi da xfoil

alpha_cl = data_cl(:, 1); % angolo di incidenza
cl_xfoil = data_cl(:, 2); % coefficiente di portanza


figure;
plot( alpha_cl, CL, 'o-', 'LineWidth',1.5);
grid on
hold on
plot(data_cl(:,1),data_cl(:,2),'LineWidth',1.5);
legend ('CL(\alpha)-matlab','CL(\alpha)-xfoil');
title ('CL(\alpha)-matlab vs CL(\alpha)-xfoil');
xlabel('\alpha', 'FontWeight','bold'); ylabel('C_L', 'FontWeight','bold');
saveas(gcf, 'Curva_cl_alpha.jpg')

