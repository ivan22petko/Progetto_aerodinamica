clc
close all
clear all

%Il presente script implementa un metodo per la stima della
%curva polare dell'ala del velivolo Cessna 172 Skyhawk, considerata come
%isolata, e ne fornisce un confronto con la stessa curva ma di un'ala con
%prestazioni ottime.
% L'ala è stata divisa in 20 settori nella direzione della corda e
%40 nella direzione dell'apertura alare. Tramite il metodo di Weissinger
%è stata calcolata, per sei angoli di incidenza tra 0 e 10 gradi, la distribuzione di 
%circolazione sull'ala, quindi, sfruttando il teorema di Kutta-Joukowsky, 
%viene calcolata la portanza complessiva dell'ala, sommando i contributi di 
%tutti i 40 settori longitudinali. 
%Per quanto rigurarda il drag indotto, è stato calcolato a partire 
%dai valori degli angoli di incidenza indotti. Questi ultimi sono stati
%calcolati per ogni settore longitudinale, tramite il calcolo della
%velocità indotta da tutti i pannelli in punti posizionati a un quarto
%della corda dell'ala. Il drag indotto complessivo è la somma dei
%contributi di tutti i 40 settori longitudinali.

%% Dati 

U_Inf_Mag = 59.2;
beta = 0;
rho = 1.225;

%% Costruzione della struttura

config.NCorpi = 1;

config.RootChord = [1.625];
config.DihedralAngle = [3]; % [°]
config.SweepAngle = [0]; % [°]
config.TaperRatio = [0.672]; 
% config.AspectRatio = [0]; 
config.Span = [11];
config.LEPosition_X = [0];
config.LEPosition_Y = [0];
config.LEPosition_Z = [0];

config.RotationAngle_X = [0];
config.RotationAngle_Y = [0];
config.RotationAngle_Z = [0];

% Discretization options
config.SemiSpanwiseDiscr = [20];
config.ChordwiseDiscr = [20];


%% Preliminary computations

% Computing the span
config.SemiSpan = config.Span./2;
% Computing the surface
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
% Computing the Tip chord
config.TipChord = config.RootChord .* config.TaperRatio;

% Compute MAC
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

%% Create the geometry structure

ControlPoints = cell(config.NCorpi, 1);
InducedPoints = cell(config.NCorpi, 1);
Normals = cell(config.NCorpi, 1);
InfiniteVortices = cell(config.NCorpi, 1);
Vortices = cell(config.NCorpi, 1);
internalMesh = cell(config.NCorpi, 1);
WingExtremes = cell(config.NCorpi, 1);

for iCorpo = 1:config.NCorpi
    [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);
end
    

%% Matrices initialization

NPanelsTot = 2 * config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
matriceA = zeros(NPanelsTot, NPanelsTot);
TermineNoto = zeros(NPanelsTot, 1);

%% Construction of the matrix

rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
   
            columnIndex = 0;
            
            ControlPointHere = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            
            
            for jCorpo = 1:config.NCorpi
                
                % Cycle on all of its chordwise panels
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                    % Cycle on all of its spanwise panels
                    for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                        
                        % Update column index
                        columnIndex = columnIndex + 1;
                        
                        % Compute the influence induced by first
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                        U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        % Compute the influence induced by finite vortex
                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        % Compute the influence induced by second
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        matriceA(rowIndex, columnIndex) = dot(U, NormalHere);
                       
                        
                    end
                end

            end
            
        
            
        end
    end
end

%% Costruzione di una matrice che contiene le coordinate dei punti nel quarto della corda di ogni settore longitudinale.
    
quarter_chord_line = cell(config.NCorpi, 1);
quarter_chord_panel = cell(config.NCorpi, 1);
    
for iCorpo = 1 : config.NCorpi
    quarter_chord_line{iCorpo} = zeros(1, 2*config.SemiSpanwiseDiscr(iCorpo)+1);
    
    %per la semi ala sinistra, viene calcolato il punto a un quarto della corda di
    %ogni linea tra due pannelli adiacenti, considerando le "tip" dei pannelli contenuti in internalMesh
    
    for j = 1 : config.SemiSpanwiseDiscr(iCorpo) 
        LE = internalMesh{iCorpo}{1, j}.LEtip(1);
        TE = internalMesh{iCorpo}{config.ChordwiseDiscr(iCorpo), j}.TEtip(1);
        quarter_chord_line{iCorpo}(1, j) = (TE-LE)/4 + LE;
    end
    
    %per la linea centrale che divide le due semiali il calcolo viene
    %effettuato separatamente, considerando la "root" dell'ultimo pannello della semi
    %ala sinistra
    
    LE = internalMesh{iCorpo}{1, j}.LERoot(1);
    TE = internalMesh{iCorpo}{config.ChordwiseDiscr(iCorpo), j}.TERoot(1);
    quarter_chord_line{iCorpo}(1, j+1) = (TE-LE)/4 + LE;

    %per la semi ala destra, viene calcolato il punto a un quarto della corda di
    %ogni linea tra due pannelli adiacenti, considerando
    %le "tip" dei pannelli contenuti in internalMesh

    for j = config.SemiSpanwiseDiscr(iCorpo)+1 : 2*config.SemiSpanwiseDiscr(iCorpo)
        LE = internalMesh{iCorpo}{1, j}.LEtip(1);
        TE = internalMesh{iCorpo}{config.ChordwiseDiscr(iCorpo), j}.TEtip(1);
        quarter_chord_line{iCorpo}(1, j+1) = (TE-LE)/4 + LE;
    end
end
          
%viene completato il calcolo dei punti a un quarto della corda per ogni
%settore longitudinale, con l'aggiunta delle altre due coordinate. 
%I punti ottenuti sono quelli in cui verrà calcolata la velocità
%indotta, ai fini del calcolo dell'angolo di incidenza indotto.

quarter_chord_panel{iCorpo} = zeros(3, 2*config.SemiSpanwiseDiscr(iCorpo));

%salvataggio delle x, y, z
for iCorpo = 1 : config.NCorpi
    for  i = 1 : 2*config.SemiSpanwiseDiscr(iCorpo)
        quarter_chord_panel{iCorpo}(1, i) = (quarter_chord_line{iCorpo}(1, i) + quarter_chord_line{iCorpo}(1, i+1))/2;
        quarter_chord_panel{iCorpo}(2, i) = (internalMesh{iCorpo}{1, i}.LEtip(2) + internalMesh{iCorpo}{1, i}.LERoot(2) )/2;
        quarter_chord_panel{iCorpo}(3, i) = (internalMesh{iCorpo}{1, i}.LEtip(3) + internalMesh{iCorpo}{1, i}.LERoot(3) )/2;
    end
end
%% Ciclo per costruzione della polare

alpha = linspace(0, 10, 6);
CL3D_polar = zeros(length(alpha), config.NCorpi);
CD3D_polar = zeros(length(alpha), config.NCorpi);

Gamma_alfa2 = cell(config.NCorpi, 1);
for iCorpo = 1:config.NCorpi
    Gamma_alfa2{iCorpo} = zeros( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );
end

%inizio del ciclo che permette di valutare più angoli di incidenza, per
%costruire la polare.

for alpha_index=1:length(alpha)

rowIndex = 0;
U_Inf = [cosd(beta)*cosd(alpha(alpha_index)) sind(beta)*cosd(alpha(alpha_index)) sind(alpha(alpha_index))] .* U_Inf_Mag;

    for iCorpo = 1:config.NCorpi
        
        % Cycle on all of its chordwise panels
        for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
            % Cycle on all of its spanwise panels
            for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
                
                % Update row index
                rowIndex = rowIndex+1;
                NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
                TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);
                
            end
        end
    end
    
    % Solve the linear system
    
    Solution = linsolve(matriceA, TermineNoto);
    
    Gamma = cell(config.NCorpi, 1);
    
    rowIndex = 0;
    for iCorpo = 1:config.NCorpi
        
        Gamma{iCorpo} = zeros( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );
        
         % Cycle on all of its chordwise panels
        for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
            % Cycle on all of its spanwise panels
            for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
                
                % Update row index
                rowIndex = rowIndex + 1;
                Gamma{iCorpo}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
                if alpha_index == 2 %salvo in un cell array indipendent i valori di gamma a incidenza alfa = 2°                    
                    Gamma_alfa2{iCorpo}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
                end

            end
        end
    end
    
    %Calcolo delle velocità indotte dai vortici semi infiniti di tutti i pannelli nei punti
    %a un quarto della corta di ciascun settore longitudinale.
    
    U_induced_tab = cell(config.NCorpi, 1);
    
    rowIndex = 0;
    
    for iCorpo = 1:config.NCorpi
        U_induced_tab{iCorpo} = zeros(3, 2*config.SemiSpanwiseDiscr(iCorpo));
        
            % Cycle on all of its spanwise panels
            for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
                
                % Update row index
                rowIndex = rowIndex + 1;
       
                columnIndex = 0;
                
                ControlPointHere = (quarter_chord_panel{iCorpo}(:, SpanPanel_i))'; %trasposizione per garantire compatibilità delle dimensioni
                NormalHere = Normals{iCorpo}{1, SpanPanel_i}.Coords; 
                
                U_induced = 0;
                
                for jCorpo = 1:config.NCorpi
                    
                    % Cycle on all of its chordwise panels
                    for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                        % Cycle on all of its spanwise panels
                        for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                            
                            % Update column index
                            columnIndex = columnIndex + 1;
                            
                            % Compute the influence induced by first
                            % semi-infinite vortex
                            Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                            Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                            U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
                                               
                            % Compute the influence induced by second
                            % semi-infinite vortex
                            Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                            Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                            U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
                            
                            U_induced = U_induced + U*Gamma{jCorpo}(ChordPanel_j, SpanPanel_j);                       
                                                 
                            
                        end
                    end
    
                end
                U_induced_tab{iCorpo}(:, SpanPanel_i) = U_induced; 
                
            end
    
    end
    
    %calcolo di 2D e 3D lift
    
    L2D_tab = cell(config.NCorpi, 1);
    L3D_tab = cell(config.NCorpi, 1);
    CL3D_tab = cell(config.NCorpi, 1);

    
    for iCorpo = 1:config.NCorpi
        delta_b = (config.Span(iCorpo)) / (config.SemiSpanwiseDiscr(iCorpo)*2);
        L2D_tab{iCorpo} = zeros(1, 2*config.SemiSpanwiseDiscr(iCorpo));
    
        for i=1:2*config.SemiSpanwiseDiscr(iCorpo)
            L2D=0;
            for j=1:config.ChordwiseDiscr(iCorpo)
                L2D= L2D + rho*U_Inf_Mag*Gamma{iCorpo}(j, i)*cosd(config.DihedralAngle(iCorpo));
            end
            L2D_tab{iCorpo}(1, i) = L2D;
        end 
        L3D_tab{iCorpo} = sum(L2D_tab{iCorpo}(1, :))*delta_b;
        CL3D_tab{iCorpo} = L3D_tab{iCorpo} / (0.5*rho*U_Inf_Mag^2*config.Surface(iCorpo));
        CL3D_polar (alpha_index,iCorpo)=  CL3D_tab{iCorpo};
    end

    %Drag
    
    %calcolo alpha indotti
    alpha_ind_tab = cell(config.NCorpi, 1);
    for iCorpo = 1 : config.NCorpi
        alpha_ind_tab{iCorpo} = zeros(1, 2*config.SemiSpanwiseDiscr(iCorpo));
    
        for i = 1 : 2*config.SemiSpanwiseDiscr(iCorpo)
            NormalHere = Normals{iCorpo}{1, i}.Coords;
            alpha_ind = -atan((dot(U_induced_tab{iCorpo}(:, i), NormalHere)/U_Inf_Mag)); %abbiamo messo il meno per far uscire il CD3D positivo
            alpha_ind_tab{iCorpo}(1, i) = alpha_ind;
        end
    end
    
    %calcolo D2D
    
    D2D_tab = cell(config.NCorpi, 1);
    D3D_tab = cell(config.NCorpi, 1);
    CD3D_tab = cell(config.NCorpi, 1);
    for iCorpo = 1 : config.NCorpi
        delta_b = (config.Span(iCorpo)) / (config.SemiSpanwiseDiscr(iCorpo)*2);
        D2D_tab{iCorpo} = zeros(1, 2*config.SemiSpanwiseDiscr(iCorpo));
        
        for i = 1 : 2*config.SemiSpanwiseDiscr(iCorpo)
            D2D_tab{iCorpo}(1, i) = L2D_tab{iCorpo}(1, i)*sin(alpha_ind_tab{iCorpo}(1, i));
        end
    
        D3D_tab{iCorpo} = sum(D2D_tab{iCorpo}(1, :))*delta_b;
        CD3D_tab{iCorpo} = D3D_tab{iCorpo} / (0.5*rho*U_Inf_Mag^2*config.Surface(iCorpo));
        CD3D_polar (alpha_index,iCorpo)= CD3D_tab{iCorpo}; 
    end

end

save('CD_cessna_polar_wing.mat', 'CD3D_polar')
save('CL_cessna_polar_wing.mat', 'CL3D_polar')

%% CL/alpha e polare 

figure(1)

plot(alpha,CL3D_polar(:,1),'-o r','LineWidth', 2);
xlabel ('\alpha [°]', 'FontWeight', 'bold');
ylabel ('C_L', 'FontWeight', 'bold');
title ('Curva Cl/alpha ala 3D');
grid on

figure(2)

plot (CD3D_polar,CL3D_polar(:,1),'-* b','LineWidth',2);
xlabel ('C_D', 'FontWeight', 'bold');
ylabel ('C_L', 'FontWeight', 'bold');
title ('Curva polare ala 3D');
grid on



%% Plot 3D dell'ala

figure(3)
hold on
colormap('jet')
colorbar

if config.NCorpi>1
    clim([min(min(min(Gamma{1})),min(min(Gamma{2}))) max(max(max(Gamma{1})), max(max(Gamma{2})))]);
else
    clim([min(min(Gamma{1})) max(max(Gamma{1}))]);
end

for iCorpo = 1:config.NCorpi
    
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            Estremo_1 = internalMesh{iCorpo}{ChordPanel_i, SpanPanel_i}.LEtip;
            Estremo_2 = internalMesh{iCorpo}{ChordPanel_i, SpanPanel_i}.TEtip;
            Estremo_3 = internalMesh{iCorpo}{ChordPanel_i, SpanPanel_i}.TERoot;
            Estremo_4 = internalMesh{iCorpo}{ChordPanel_i, SpanPanel_i}.LERoot;
            
            X = [Estremo_1(1), Estremo_2(1), Estremo_3(1), Estremo_4(1)];
            Y = [Estremo_1(2), Estremo_2(2), Estremo_3(2), Estremo_4(2)];
            Z = [Estremo_1(3), Estremo_2(3), Estremo_3(3), Estremo_4(3)];
            
            
            patch('XData',X,'YData', Y,'ZData', Z, "FaceColor", 'flat','FaceVertexCData', Gamma{iCorpo}(ChordPanel_i,SpanPanel_i), 'EdgeColor', 'k');
            
            
        end
    end
end

axis equal
view(3)
grid minor
xlabel('X', 'FontWeight', 'bold');
ylabel('Y', 'FontWeight', 'bold');
zlabel('Z', 'FontWeight', 'bold', 'Rotation', 0);
title('Distribuzione di Gamma sulla superficie');
set(colorbar, 'position', [0.85,0.1,0.03,0.8])

%% PLOT CIRCOLAZIONE per alpha=2
%nel grafico viene rappresentata la distribuzione di circolazione nell'ultima
%striscia trasversale.

z=linspace (-config.SemiSpan(1),config.SemiSpan(1),2*config.SemiSpanwiseDiscr(1));
gamma_plot_cessna = Gamma_alfa2{1}(1, :);
save('gamma_plot_cessna.mat', 'gamma_plot_cessna')

figure(4)

plot(z,gamma_plot_cessna, 'b-', 'LineWidth',2);
grid on
xlabel('z', 'FontWeight', 'bold')
ylabel('\Gamma', 'FontWeight', 'bold')
title('Distribuzione di gamma sull''apertura')

pbaspect([2.5 1 1]); 

%% Confronto con ala ellittica

% Dichiarazione dati ala ellittica, imponiamo che il coefficiente di
% portanza prodotto dall'ala ellittica sia lo stesso di quello ottenuto dal
% metodo implementato, e come parametri geometrici chiediamo la stessa
% superficie e la stessa apertura alare. 

cl_alpha = 2*pi; % ipotizzati
alpha_0 = 0; % ipotizzata

U_inf = U_Inf_Mag;
S = config.Surface;
b = config.Span; 
c0 = (4*S)/(pi*b); 
lambda = b^2/S; %aspect ratio ala ellittica

CL_alpha = (cl_alpha) / (1+(cl_alpha)/pi*lambda);  %calcolo della pendenza della curva CL-alpha dell'ala ellittica


CD_ell = [];

for i = 1 : length(alpha)
    CD = (CL3D_polar(i)^2)/(pi*lambda);

    CD_ell = [CD_ell; CD];

end

%% plot dei confronti

%plot polare dell'ala ellittica
figure
plot (CD_ell, CL3D_polar,'-* r','LineWidth',1);
hold on
plot (CD3D_polar, CL3D_polar,'-o b','LineWidth',1);
xlabel ('C_D', 'FontWeight', 'bold');
ylabel ('C_L', 'FontWeight', 'bold');
title ('Confronto polare ala ellittica e cessna');
legend ('Ellittica', 'Cessna', 'Location', 'best')
grid on

%calcolo della efficienza 
E_cessna = CL3D_polar(2)/CD3D_polar(2);
E_ellittica = CL3D_polar(2)/CD_ell(2); 

%% confronto andamenti circolazione con alpha 2 gradi tra ala ellittica e ala Cessna
alpha = deg2rad(2); 

%ala ellittica
gamma0 = ((2*U_inf*pi*c0*b) / (2*b+ pi*c0))*(deg2rad(alpha_0)-alpha);
z=linspace(-b/2, b/2, 40);
gamma_plot_ell = gamma0 * sqrt(1 - (z.^2)/(b/2)^2);
save('gamma_plot_ell.mat', 'gamma_plot_ell')

%cessna
z_cessna=linspace (-config.SemiSpan(1),config.SemiSpan(1),2*config.SemiSpanwiseDiscr(1));
gamma_plot_cessna = Gamma_alfa2{1}(1, :);

figure
plot(z, abs(gamma_plot_ell), '-* r','LineWidth',1)
hold on
plot(z_cessna, abs(gamma_plot_cessna), '-o b','LineWidth',1)
grid on
xlabel('z', 'FontWeight', 'bold')
ylabel('\Gamma', 'FontWeight', 'bold', 'Rotation', 0)
title('Confronto andamenti circolazione')
legend('Ala ellittica', 'Cessna', 'Location', 'best')
pbaspect([2.75 1 1]);

