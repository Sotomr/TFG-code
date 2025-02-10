clc;
clear;
close all;

% Rangs de paràmetres
porcentajeOpt_values = 0:10:100;        % Percentatge d'optimització
porcentajeGroup_values = 0:10:100;       % Percentatge de grups
maxNoImprovement_values = 5:5:10;       % Nombre màxim d'iteracions sense millora

% Paràmetres fixos per a gràfiques
porcentajeGroup_fixed = 80;
porcentajeOpt_fixed_values = [0, 25, 50, 100];
numIteracions_fixed = 300;
maxNoImprovement_fixed = 1;

%% 1. Gràfica: Cost i Satisfacció vs Percentatge d'Optimització
disp('Generant gràfiques de Cost i Satisfacció vs Percentatge d''Optimització...');
cost_opt = zeros(length(porcentajeOpt_values),1);
satisf_opt = zeros(length(porcentajeOpt_values),1);
for i = 1:length(porcentajeOpt_values)
    opt = porcentajeOpt_values(i);
    grp = porcentajeGroup_fixed;
    [cost, satisf, ~] = EQMOz(opt, grp, numIteracions_fixed, maxNoImprovement_fixed,100,25);
    cost_opt(i) = cost;
    satisf_opt(i) = satisf;
end

figure('Color','w');
subplot(2,1,1);
plot(porcentajeOpt_values, cost_opt, 'r-o','LineWidth',1.5);
xlabel('Percentatge d''Optimització','FontSize',10);
ylabel('Cost Total General','FontSize',10);
title('Cost vs Percentatge d''Optimització','FontSize',10);
grid on;

subplot(2,1,2);
plot(porcentajeOpt_values, satisf_opt, 'b--o','LineWidth',1.5);
xlabel('Percentatge d''Optimització','FontSize',10);
ylabel('Satisfacció Mitjana Global','FontSize',10);
title('Satisfacció vs Percentatge d''Optimització','FontSize',10);
grid on;

%% 2. Gràfica: Cost i Satisfacció vs Percentatge de Grup
disp('Generant gràfiques de Cost i Satisfacció vs Percentatge de Grup...');
numOpt = length(porcentajeOpt_fixed_values);
cost_group = zeros(numOpt, length(porcentajeGroup_values));
satisf_group = zeros(numOpt, length(porcentajeGroup_values));
for i = 1:numOpt
    opt_fixed = porcentajeOpt_fixed_values(i);
    for j = 1:length(porcentajeGroup_values)
        grp = porcentajeGroup_values(j);
        [cost, satisf, ~] = EQMOz(opt_fixed, grp, numIteracions_fixed, maxNoImprovement_fixed,100,25);
        cost_group(i,j) = cost;
        satisf_group(i,j) = satisf;
    end
end

figure('Color','w');
for i = 1:numOpt
    subplot(2,1,1);
    plot(porcentajeGroup_values, cost_group(i,:), 'LineWidth',1.5);
    hold on;
    subplot(2,1,2);
    plot(porcentajeGroup_values, satisf_group(i,:), 'LineWidth',1.5);
    hold on;
end
subplot(2,1,1);
xlabel('Percentatge de Grup','FontSize',10);
ylabel('Cost Total General','FontSize',10);
title('Cost vs Percentatge de Grup','FontSize',10);
legend(arrayfun(@(x) sprintf('%d%% Opt.', x), porcentajeOpt_fixed_values, 'UniformOutput', false));
grid on;
hold off;
subplot(2,1,2);
xlabel('Percentatge de Grup','FontSize',10);
ylabel('Satisfacció Mitjana Global','FontSize',10);
title('Satisfacció vs Percentatge de Grup','FontSize',10);
legend(arrayfun(@(x) sprintf('%d%% Opt.', x), porcentajeOpt_fixed_values, 'UniformOutput', false));
grid on;
hold off;

%% 3. Matrius de Correlació: una per Satisfacció i una per Cost
disp('Generant matrius de correlació...');

corr_satisf = zeros(length(porcentajeOpt_values), length(porcentajeGroup_values));
corr_cost   = zeros(length(porcentajeOpt_values), length(porcentajeGroup_values));

for i = 1:length(porcentajeOpt_values)
    for j = 1:length(porcentajeGroup_values)
        opt = porcentajeOpt_values(i);
        grp = porcentajeGroup_values(j);
        [cost, satisf, ~] = EQMOz(opt, grp, numIteracions_fixed, maxNoImprovement_fixed,100,25);
        corr_satisf(i,j) = satisf;
        corr_cost(i,j)   = cost;
    end
end


figure('Color','w');
imagesc(porcentajeGroup_values, porcentajeOpt_values, corr_satisf);
colorbar;
xlabel('% de Grup','FontSize',10);
ylabel('% d''Optimització','FontSize',10);
title('Matriu de Correlació: Satisfacció','FontSize',10);
colormap(gca, [linspace(0,0.2,256)', linspace(0,1,256)', linspace(0,0.2,256)']);
set(gca, 'YDir', 'normal');

figure('Color','w');
imagesc(porcentajeGroup_values, porcentajeOpt_values, corr_cost);
colorbar;
xlabel('% de Grup','FontSize',10);
ylabel('% d''Optimització','FontSize',10);
title('Matriu de Correlació: Cost','FontSize',10);
colormap(gca, [linspace(0,1,256)', zeros(256,1), zeros(256,1)]);
set(gca, 'YDir', 'normal');

%% 4. Gràfica: Cost i Satisfacció vs Nombre Màxim d'Iteracions sense Millora
disp('Generant gràfiques de Cost i Satisfacció vs Nombre Màxim d''Iteracions sense Millora...');
cost_iterations = zeros(length(maxNoImprovement_values),1);
satisf_iterations = zeros(length(maxNoImprovement_values),1);
for i = 1:length(maxNoImprovement_values)
    maxNoImp = maxNoImprovement_values(i);
    opt = porcentajeOpt_fixed_values(end);  % 100%
    grp = porcentajeGroup_fixed;
    [cost, satisf, ~] = EQMOz(opt, grp, numIteracions_fixed, maxNoImp,100,25);
    cost_iterations(i) = cost;
    satisf_iterations(i) = satisf;
end

figure('Color','w');
yyaxis left
plot(maxNoImprovement_values, cost_iterations, 'r-o','LineWidth',1.5);
ylabel('Cost Total General','FontSize',10);
yyaxis right
plot(maxNoImprovement_values, satisf_iterations, 'b--o','LineWidth',1.5);
ylabel('Satisfacció Mitjana Global','FontSize',10);
xlabel('Nombre Màxim d''Iteracions sense Millora','FontSize',10);
title('Cost i Satisfacció vs Iteracions sense Millora','FontSize',10);
legend({'Cost Total General','Satisfacció Mitjana Global'},'Location','northwest');
grid on;
