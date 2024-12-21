clear all
close all
clc

%Il presente codice sfrutta i risultati dei codici che calcolano la polare
%di Cessna 172 e BAE 100 Hawk per calcolare la polare dei due velivoli
%completi ed effettuare un confronto grafico tra la polare del velivolo
%completo e del velivolo con ala isolata.

%% CESSNA

load('CD_cessna_polar_wing.mat', 'CD3D_polar');
CD_cessna_polar_wing = CD3D_polar;
load('CL_cessna_polar_wing.mat', 'CL3D_polar');
CL_cessna_polar_wing = CL3D_polar;

load('CD_cessna_polar_wingtail.mat', 'CD_cessna_polar_wingtail');
load('CL_cessna_polar_wingtail.mat', 'CL_cessna_polar_wingtail');

figure(1)
plot(CD_cessna_polar_wing, CL_cessna_polar_wing, '--* r', 'LineWidth', 1)
hold on
plot(CD_cessna_polar_wingtail, CL_cessna_polar_wingtail, '-o r', 'LineWidth', 1)
grid on
% legend('Polare ala Cessna 172', 'Polare ala+coda Cessna 172', 'Location', 'southeast')
% title('Confronto polari Cessna 172')
xlabel('C_D', 'FontWeight', 'bold'); ylabel('C_L', 'FontWeight', 'bold')

%% HAWK 

load('CD_hawk_polar_wing.mat', 'CD3D_polar');
CD_hawk_polar_wing = CD3D_polar;
load('CL_hawk_polar_wing.mat', 'CL3D_polar');
CL_hawk_polar_wing = CL3D_polar;

load('CD_hawk_polar_wingtail.mat', 'CD_hawk_polar_wingtail');
load('CL_hawk_polar_wingtail.mat', 'CL_hawk_polar_wingtail');

% figure(2)
plot(CD_hawk_polar_wing, CL_hawk_polar_wing, '--* b', 'LineWidth', 1')
hold on
plot(CD_hawk_polar_wingtail, CL_hawk_polar_wingtail, '-o b', 'LineWidth', 1)
grid on
legend('Polare ala Cessna 172', 'Polare ala+coda Cessna 172','Polare ala Hawk 100', 'Polare ala+coda Hawk 100', 'Location', 'southeast')
% title('Confronto polari')
xlabel('C_D', 'FontWeight', 'bold'); ylabel('C_L', 'FontWeight', 'bold')

saveas(gcf, 'confronto_polari.jpg')