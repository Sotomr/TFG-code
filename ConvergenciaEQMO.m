%% Script per analitzar el temps de convergència de EQMO2 i sobreposar tendències teòriques

%   O(N): Creixement lineal
%   O(N log N): Creixement quasi-lineal
%   O(N^2): Creixement quadràtic
%


clc;
clear;
close all;

%% Paràmetres per a EQMO2
porcentajeOpt    = 50;       % Percentatge d'optimització (exemple)
porcentajeGroup  = 50;       % Percentatge de grup (exemple)
numIteraciones   = 300;      % Nombre màxim d'iteracions
maxNoImprovement = 10;       % Nombre màxim d'iteracions sense millora

%% Definició de les situacions
numScenarios = 20;  
numP_array = round(linspace(10, 250, numScenarios));  % de 10 a 170 paquets
numD_array = round(linspace(2, 40, numScenarios));     % de 2 a 40 drons

labels = cell(1, numScenarios);        % Etiquetes de cada 
time_values = zeros(1, numScenarios);  % Temps de convergència (segons)

%% Mesurar el temps de convergència per cada situació
for i = 1:numScenarios
    numP = numP_array(i);
    numD = numD_array(i);
    labels{i} = sprintf('%dP-%dD', numP, numD);
    fprintf('Executant EQMO2 per la situació %s...\n', labels{i});
    tic;
    [~, ~, bestIter] = EQMO(porcentajeOpt, porcentajeGroup, numIteraciones, maxNoImprovement, numP, numD);
    time_values(i) = toc;
    fprintf('Situació %s: Temps de convergència = %.4f segons.\n', labels{i}, time_values(i));
end

%% Ajustar models de tendència sobre la base de les dades experimentals
x = 1:numScenarios;  % Índex numèric per a cada situació

% Model 1: O(N) - ajust lineal
coeff_lin = polyfit(x, time_values, 1);
y_lin = polyval(coeff_lin, x);

% Model 2: O(N log N) - ajust lineal sobre x*log(x)
X_log = [ones(numScenarios,1), (x'.*log(x'))];
coeff_nlogn = X_log \ time_values';
y_nlogn = coeff_nlogn(1) + coeff_nlogn(2)*(x.*log(x));

% Model 3: O(N^2) - ajust quadràtic
coeff_quad = polyfit(x, time_values, 2);
y_quad = polyval(coeff_quad, x);

%% Crear la gràfica amb els ajustaments sol·licitats
figure('Color','w','Position',[100, 100, 900, 600]);

p_exp = plot(x, time_values, 'r-', 'LineWidth', 2, ...
    'Marker','o', 'MarkerSize',4, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
hold on;

p_lin = plot(x, y_lin, 'b--', 'LineWidth', 2);
p_nlogn = plot(x, y_nlogn, 'g--', 'LineWidth', 2);
p_quad = plot(x, y_quad, 'm--', 'LineWidth', 2);

grid on;
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', x, 'XTickLabel', labels, 'FontSize',9);
xlabel('Situació (Paquets - Drons)', 'FontSize',11, 'FontWeight','bold');
ylabel('Temps de Convergència (segons)', 'FontSize',11, 'FontWeight','bold');
title('Convergència de EQMO i Tendències de Creixement', 'FontSize',14, 'FontWeight','bold');
xlim([min(x)-0.5, max(x)+0.5]);
ylim([min(time_values)*0.9, max(time_values)*1.1]);
box on;

legend([p_exp, p_lin, p_nlogn, p_quad], ...
    {'Dades experimentals', 'O(N): Creixement lineal', 'O(N log N): Creixement quasi-lineal', 'O(N^2): Creixement quadràtic'}, ...
    'Location', 'northwest', 'FontSize',10);
