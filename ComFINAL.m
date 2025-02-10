%% Script de Gràfiques per comparar els 5 models EQMO
clc;
clear;
close all;

% Rangs de paràmetres
porcentajeOpt_values = 0:10:100;        % Percentatge d'optimització
porcentajeGroup_fixed = 80;             % Percentatge de grups fix
numIteracions_fixed = 300;              % Número d'iteracions fix
maxNoImprovement_fixed = 1;             % Nombre màxim d'iteracions sense millora

nOpt = length(porcentajeOpt_values);


cost_EQMO_base   = zeros(nOpt,1);
satisf_EQMO_base = zeros(nOpt,1);

cost_EQMO    = zeros(nOpt,1);      
satisf_EQMO  = zeros(nOpt,1);

cost_EQMOv1   = zeros(nOpt,1);
satisf_EQMOv1 = zeros(nOpt,1);

cost_EQMOv2   = zeros(nOpt,1);
satisf_EQMOv2 = zeros(nOpt,1);

cost_EQMOv3   = zeros(nOpt,1);
satisf_EQMOv3 = zeros(nOpt,1);

disp('Generant gràfiques de Cost i Satisfacció vs Percentatge d''Optimització per als 5 models...');

for i = 1:nOpt
    opt = porcentajeOpt_values(i);
    grp = porcentajeGroup_fixed;
    
    % 5 models
    [cost0, satisf0, ~] = EQMO(opt, grp, numIteracions_fixed, maxNoImprovement_fixed);
    [cost1, satisf1, ~] = EQMOy(opt, grp, numIteracions_fixed, maxNoImprovement_fixed,100,25);
    [cost2, satisf2, ~] = EQMOv1(opt, grp, numIteracions_fixed, maxNoImprovement_fixed,100,25);
    [cost3, satisf3, ~] = EQMOv2(opt, grp, numIteracions_fixed, maxNoImprovement_fixed,100,25);
    [cost4, satisf4, ~] = EQMOv3(opt, grp, numIteracions_fixed, maxNoImprovement_fixed,100,25);
    
    cost_EQMO_base(i)   = cost0;
    satisf_EQMO_base(i) = satisf0;
    
    cost_EQMO(i)    = cost1;
    satisf_EQMO(i)  = satisf1;
    
    cost_EQMOv1(i)   = cost2;
    satisf_EQMOv1(i) = satisf2;
    
    cost_EQMOv2(i)   = cost3;
    satisf_EQMOv2(i) = satisf3;
    
    cost_EQMOv3(i)   = cost4;
    satisf_EQMOv3(i) = satisf4;
end



figure('Color','w');


subplot(2,1,1);
plot(porcentajeOpt_values, cost_EQMO_base, 'm-o','LineWidth',1.5); hold on;
plot(porcentajeOpt_values, cost_EQMO, 'r-o','LineWidth',1.5);
plot(porcentajeOpt_values, cost_EQMOv1, 'b-s','LineWidth',1.5);
plot(porcentajeOpt_values, cost_EQMOv2, 'g-^','LineWidth',1.5);
plot(porcentajeOpt_values, cost_EQMOv3, 'k--d','LineWidth',1.5);
xlabel('Percentatge d''Optimització','FontSize',10);
ylabel('Cost Total General','FontSize',10);
title('Cost vs Percentatge d''Optimització','FontSize',10);
legend('EQMO','EQMOy','EQMOv1','EQMOv2','EQMOv3','Location','best');
grid on;

subplot(2,1,2);
plot(porcentajeOpt_values, satisf_EQMO_base, 'm-o','LineWidth',1.5); hold on;
plot(porcentajeOpt_values, satisf_EQMO, 'r-o','LineWidth',1.5);
plot(porcentajeOpt_values, satisf_EQMOv1, 'b-s','LineWidth',1.5);
plot(porcentajeOpt_values, satisf_EQMOv2, 'g-^','LineWidth',1.5);
plot(porcentajeOpt_values, satisf_EQMOv3, 'k--d','LineWidth',1.5);
xlabel('Percentatge d''Optimització','FontSize',10);
ylabel('Satisfacció Mitjana Global','FontSize',10);
title('Satisfacció vs Percentatge d''Optimització','FontSize',10);
legend('EQMO','EQMOy','EQMOv1','EQMOv2','EQMOv3','Location','best');
grid on;
