clear
close all
clc

%Il seguente codice permette di generare un grafico con l'andamento ella
%circolazione (ad incidenza 2 gradi) dell'ala del Cessna 172, del BAE
%Systems Hawk 100 (considerate come isolate) e dell'ala ellittica
%utilizzata per confrontare le prestazioni dell'ala del Cessna con quelle
%di un'ala ottima.

%%
config.SemiSpanwiseDiscr = [20];
config.ChordwiseDiscr = [20];

config.Span = [11];
config.SemiSpan = config.Span./2;
z_cessna=linspace (-config.SemiSpan(1),config.SemiSpan(1),2*config.SemiSpanwiseDiscr(1));
load('gamma_plot_cessna.mat', 'gamma_plot_cessna');

b = config.Span;
z_ell=linspace(-b/2, b/2, 40);
load('gamma_plot_ell.mat', 'gamma_plot_ell');

config.Span = [9.08];
config.SemiSpan = config.Span./2;
z_hawk=linspace (-config.SemiSpan(1),config.SemiSpan(1),2*config.SemiSpanwiseDiscr(1));
load('gamma_plot_hawk.mat', 'gamma_plot_hawk');

figure
plot(z_cessna, gamma_plot_cessna, 'LineWidth', 2)
hold on
plot(z_ell, abs(gamma_plot_ell), 'LineWidth', 2)
plot(z_hawk, gamma_plot_hawk, 'LineWidth', 2)
grid on
xlabel('Z', 'FontWeight', 'bold'); ylabel('\Gamma', 'Rotation', 0, 'FontWeight', 'bold')
legend('Cessna 172', 'Ala ellittica', 'Hawk 100')
pbaspect([2.75 1 1]);

%%
saveas(gcf, 'confronto_cicrcolazioni.jpg')
