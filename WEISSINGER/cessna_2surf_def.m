clc
close all
clear all

%Il presente script implementa un metodo per la stima della
%curva polare dell'ala del velivolo Cessna 172 Skyhawk, considerata come
%accoppiata al suo piano di coda.
%L'ala è stata divisa in 20 settori nella direzione della corda e
%40 nella direzione dell'apertura alare, mentre la coda in 10 e 
%20 rispettivamente. Tramite il metodo di Weissinger è stata calcolata, per sei 
%angoli di incidenza tra 0 e 10 gradi, la distribuzione di 
%circolazione sull'ala, quindi, sfruttando il teorema di Kutta-Joukowsky, 
%viene calcolata la portanza complessiva dell'ala, sommando i contributi di 
%tutti i 40 settori longitudinali. 
%Per quanto rigurarda il drag indotto, è stato calcolato a partire 
%dai valori degli angoli di incidenza indotti. Questi ultimi sono stati
%calcolati per ogni settore longitudinale, tramite il calcolo della
%velocità indotta da tutti i pannelli in punti posizionati a un quarto
%della corda dell'ala. Il drag indotto complessivo è la somma dei
%contributi di tutti i 40 settori longitudinali.

%% Dati da slide

U_Inf_Mag = 59.2;
beta = 0;
rho = 1.225;

%% Costruzione della struttura

config.NCorpi = 2;

config.RootChord = [1.625, 1.4 ];
config.DihedralAngle = [3, 0]; % [°]
config.SweepAngle = [0, 0]; % [°]
config.TaperRatio = [0.672, 0.5714]; 
% config.AspectRatio = [0]; 
config.Span = [11, 3.4];
config.LEPosition_X = [0, 4.668];
config.LEPosition_Y = [0, 0];
config.LEPosition_Z = [0, -0.2];

config.RotationAngle_X = [0, 0];
config.RotationAngle_Y = [0, 0];
config.RotationAngle_Z = [0, 0];

% Discretization options
config.SemiSpanwiseDiscr = [20, 10];
config.ChordwiseDiscr = [20, 10];

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

NPanelsTot = 2* config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
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

alpha = linspace(0,10,6);
L3D_wingtail = zeros(6, 2);
D3D_wingtail = zeros(6, 2);
CD_cessna_polar_wingtail = zeros(6, 1);
CL_cessna_polar_wingtail = zeros(6, 1);
check = 1;
CL3D_polar = zeros(length(alpha),config.NCorpi);
CD3D_polar = zeros(length(alpha),config.NCorpi);

Gamma_alfa2 = cell(config.NCorpi, 1);
for iCorpo = 1:config.NCorpi
    Gamma_alfa2{iCorpo} = zeros( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );
end

for alpha_index=1:length (alpha)

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
                
                ControlPointHere = (quarter_chord_panel{iCorpo}(:, SpanPanel_i) )'; %traspongo per compatibilità delle dimensioni
                NormalHere = Normals{iCorpo}{1, SpanPanel_i}.Coords; %rimane la stessa del pannello che è quella che aveva scritto il prof
                
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
                U_induced_tab{iCorpo}(:, SpanPanel_i) = U_induced; % u nei punti di quarto di corda 
                
            end
    
    end
    
    %calcolo di 2D e 3D Lift
    
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
        CL3D_polar (alpha_index,iCorpo) =  CL3D_tab{iCorpo};
        L3D_wingtail(check, iCorpo) = L3D_tab{iCorpo};
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
        D3D_wingtail(check, iCorpo) = D3D_tab{iCorpo};
    end
    
    CD_cessna_polar_wingtail(check) = (D3D_wingtail(check, 1) + D3D_wingtail(check, 2))/(0.5*rho*U_Inf_Mag^2*config.Surface(1));
    CL_cessna_polar_wingtail(check) = (L3D_wingtail(check, 1) + L3D_wingtail(check, 2))/(0.5*rho*U_Inf_Mag^2*config.Surface(1));

    check = check + 1;

end

save('CD_cessna_polar_wingtail.mat', 'CD_cessna_polar_wingtail')
save('CL_cessna_polar_wingtail.mat', 'CL_cessna_polar_wingtail')
%% CL/alpha e polare 

figure(1)

plot(alpha,CL3D_polar(:,1),'-o r','LineWidth', 2);
xlabel ('\alpha [°]', 'FontWeight', 'bold');
ylabel ('C_L', 'FontWeight', 'bold');
title ('Curva C_L(\alpha) ala 3D');
grid on

figure(2)

plot (CD3D_polar(:,1),CL3D_polar(:,1),'-* b','LineWidth',2);
xlabel ('C_D', 'FontWeight', 'bold');
ylabel ('C_L', 'FontWeight', 'bold');
title ('Curva polare ala 3D');
grid on



%% Plot 3D di ala+coda

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
gamma_plot= Gamma_alfa2{1}(1, :);

figure(4)
grid on
plot(z,gamma_plot, '-b', 'LineWidth',2);
xlabel('z', 'FontWeight', 'bold');
ylabel('\Gamma', 'FontWeight', 'bold', 'Rotation', 0);
grid on
title('Distribuzione della circolazione sull''apertura');

pbaspect([2.5 1 1]); 