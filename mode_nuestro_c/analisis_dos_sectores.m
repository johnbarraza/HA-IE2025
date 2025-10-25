% ANALISIS_DOS_SECTORES.m
% Script para analizar resultados del modelo Aiyagari con dos sectores
% y realizar experimentos de política

clear all; clc; close all;

%% =========================================================================
% EXPERIMENTO 1: EFECTOS DE LA TASA IMPOSITIVA
% =========================================================================

fprintf('EXPERIMENTO 1: Efectos de cambios en la tasa impositiva\n');
fprintf('=========================================================\n\n');

% Tasas impositivas a evaluar
tau_vec = [0.1, 0.2, 0.3, 0.4];
n_tau = length(tau_vec);

% Almacenar resultados
results_tau = struct();

for t = 1:n_tau
    fprintf('\nEvaluando tau = %.2f\n', tau_vec(t));
    
    % Modificar el parámetro tau en el modelo principal
    % (Aquí deberías correr el modelo con cada tau)
    
    % Por ahora, simulamos algunos resultados esperados
    % En la práctica, ejecutarías aiyagari_dos_sectores.m con cada tau
    
    % Resultados simulados (reemplazar con valores reales)
    results_tau(t).tau = tau_vec(t);
    results_tau(t).share_informal_empleo = 0.20 + 0.3*tau_vec(t); % Aumenta con tau
    results_tau(t).share_informal_producto = 0.15 + 0.25*tau_vec(t);
    results_tau(t).gini = 0.35 + 0.1*tau_vec(t);
    results_tau(t).r_eq = 0.04 - 0.01*tau_vec(t);
end

% Graficar resultados
figure('Name', 'Efectos de la Tasa Impositiva', 'Position', [100 100 800 600])

subplot(2,2,1)
plot([results_tau.tau], [results_tau.share_informal_empleo], 'b-o', 'LineWidth', 2)
xlabel('Tasa Impositiva τ')
ylabel('Share Informal (%)')
title('Informalidad en Empleo')
grid on

subplot(2,2,2)
plot([results_tau.tau], [results_tau.share_informal_producto], 'r-o', 'LineWidth', 2)
xlabel('Tasa Impositiva τ')
ylabel('Share Informal (%)')
title('Informalidad en Producto')
grid on

subplot(2,2,3)
plot([results_tau.tau], [results_tau.gini], 'g-o', 'LineWidth', 2)
xlabel('Tasa Impositiva τ')
ylabel('Índice de Gini')
title('Desigualdad de Riqueza')
grid on

subplot(2,2,4)
plot([results_tau.tau], [results_tau.r_eq], 'm-o', 'LineWidth', 2)
xlabel('Tasa Impositiva τ')
ylabel('Tasa de Interés')
title('Tasa de Interés de Equilibrio')
grid on

%% =========================================================================
% EXPERIMENTO 2: DECOMPOSICIÓN DE LA INFORMALIDAD POR RIQUEZA
% =========================================================================

fprintf('\n\nEXPERIMENTO 2: Informalidad por niveles de riqueza\n');
fprintf('====================================================\n\n');

% Simulación de datos (reemplazar con resultados del modelo)
n_percentiles = 20;
percentiles = linspace(5, 100, n_percentiles);

% Patrón típico: informalidad decrece con riqueza
share_informal_by_wealth = 0.8 * exp(-0.05 * percentiles) + 0.1;
hours_total_by_wealth = 0.6 + 0.2 * exp(-0.03 * percentiles);

figure('Name', 'Informalidad por Riqueza', 'Position', [100 100 800 400])

subplot(1,2,1)
plot(percentiles, share_informal_by_wealth, 'b-', 'LineWidth', 2)
xlabel('Percentil de Riqueza')
ylabel('Share de Horas en Informal')
title('Informalidad vs Riqueza')
grid on
ylim([0 1])

subplot(1,2,2)
plot(percentiles, hours_total_by_wealth, 'r-', 'LineWidth', 2)
xlabel('Percentil de Riqueza')
ylabel('Horas Totales Trabajadas')
title('Oferta Laboral vs Riqueza')
grid on
ylim([0 1])

%% =========================================================================
% EXPERIMENTO 3: TRANSICIONES ENTRE SECTORES
% =========================================================================

fprintf('\n\nEXPERIMENTO 3: Análisis de transiciones\n');
fprintf('========================================\n\n');

% Matriz de transición simulada (4 estados: 2 productividades x 2 sectores principales)
% Estados: [z1-formal, z1-informal, z2-formal, z2-informal]

% Para agentes con baja riqueza
P_low_wealth = [
    0.70, 0.20, 0.08, 0.02;  % z1-formal
    0.15, 0.75, 0.02, 0.08;  % z1-informal
    0.25, 0.05, 0.65, 0.05;  % z2-formal
    0.02, 0.28, 0.10, 0.60   % z2-informal
];

% Para agentes con alta riqueza
P_high_wealth = [
    0.85, 0.05, 0.09, 0.01;  % z1-formal
    0.30, 0.50, 0.10, 0.10;  % z1-informal
    0.10, 0.02, 0.85, 0.03;  % z2-formal
    0.05, 0.10, 0.25, 0.60   % z2-informal
];

% Visualizar matrices de transición
figure('Name', 'Matrices de Transición', 'Position', [100 100 800 400])

subplot(1,2,1)
imagesc(P_low_wealth)
colorbar
title('Transiciones - Baja Riqueza')
xlabel('Estado Destino')
ylabel('Estado Origen')
set(gca, 'XTick', 1:4, 'XTickLabel', {'z1-F', 'z1-I', 'z2-F', 'z2-I'})
set(gca, 'YTick', 1:4, 'YTickLabel', {'z1-F', 'z1-I', 'z2-F', 'z2-I'})

subplot(1,2,2)
imagesc(P_high_wealth)
colorbar
title('Transiciones - Alta Riqueza')
xlabel('Estado Destino')
ylabel('Estado Origen')
set(gca, 'XTick', 1:4, 'XTickLabel', {'z1-F', 'z1-I', 'z2-F', 'z2-I'})
set(gca, 'YTick', 1:4, 'YTickLabel', {'z1-F', 'z1-I', 'z2-F', 'z2-I'})

%% =========================================================================
% CALIBRACIÓN SUGERIDA PARA DIFERENTES PAÍSES
% =========================================================================

fprintf('\n\nCALIBRACIÓN SUGERIDA PARA DIFERENTES ECONOMÍAS\n');
fprintf('================================================\n\n');

% Definir calibraciones para diferentes tipos de economías
calibrations = struct();

% Economía Desarrollada (ej. USA, Europa)
calibrations(1).name = 'Economía Desarrollada';
calibrations(1).tau = 0.30;      % Alta tasa impositiva
calibrations(1).psi_i = 1.5;     % Alta desutilidad informal (más enforcement)
calibrations(1).Ai = 0.05;       % Baja productividad informal
calibrations(1).al_i = 0.6;      % Retornos muy decrecientes
calibrations(1).target_informal = 0.10; % 10% informal

% Economía Emergente (ej. México, Brasil)
calibrations(2).name = 'Economía Emergente';
calibrations(2).tau = 0.20;      % Tasa impositiva media
calibrations(2).psi_i = 1.0;     % Desutilidad media
calibrations(2).Ai = 0.08;       % Productividad informal media
calibrations(2).al_i = 0.75;     % Retornos moderadamente decrecientes
calibrations(2).target_informal = 0.30; % 30% informal

% Economía en Desarrollo (ej. Colombia, Perú)
calibrations(3).name = 'Economía en Desarrollo';
calibrations(3).tau = 0.15;      % Baja tasa impositiva efectiva
calibrations(3).psi_i = 0.8;     % Baja desutilidad (poco enforcement)
calibrations(3).Ai = 0.10;       % Productividad informal relativamente alta
calibrations(3).al_i = 0.85;     % Retornos menos decrecientes
calibrations(3).target_informal = 0.45; % 45% informal

% Mostrar tabla de calibraciones
fprintf('%-25s | τ    | ψⁱ   | Aⁱ   | αⁱ   | Target Informal\n', 'Economía');
fprintf('--------------------------------------------------------\n');
for i = 1:length(calibrations)
    fprintf('%-25s | %.2f | %.2f | %.2f | %.2f | %.0f%%\n', ...
        calibrations(i).name, ...
        calibrations(i).tau, ...
        calibrations(i).psi_i, ...
        calibrations(i).Ai, ...
        calibrations(i).al_i, ...
        calibrations(i).target_informal*100);
end

%% =========================================================================
% FUNCIONES DE IMPULSO-RESPUESTA
% =========================================================================

fprintf('\n\nANÁLISIS DE IMPULSO-RESPUESTA\n');
fprintf('==============================\n\n');

% Simular respuesta a un shock de productividad en el sector informal
T = 20; % Períodos
t = 1:T;

% Shock: aumento de 10% en productividad informal en t=5
shock_period = 5;
Ai_base = 0.08;
Ai_shock = Ai_base * 1.10;

% Respuestas simuladas (en la práctica, resolver el modelo en cada período)
Li_response = zeros(T,1);
Lf_response = zeros(T,1);
r_response = zeros(T,1);

for period = 1:T
    if period < shock_period
        Li_response(period) = 0.30;
        Lf_response(period) = 0.70;
        r_response(period) = 0.040;
    elseif period == shock_period
        % Impacto inicial
        Li_response(period) = 0.35;
        Lf_response(period) = 0.65;
        r_response(period) = 0.038;
    else
        % Convergencia gradual
        decay = 0.8^(period - shock_period);
        Li_response(period) = 0.30 + 0.05*decay;
        Lf_response(period) = 0.70 - 0.05*decay;
        r_response(period) = 0.040 - 0.002*decay;
    end
end

figure('Name', 'Impulso-Respuesta', 'Position', [100 100 800 600])

subplot(3,1,1)
plot(t, Li_response, 'b-', 'LineWidth', 2)
hold on
plot([shock_period shock_period], [min(Li_response) max(Li_response)], 'r--')
xlabel('Período')
ylabel('L^i')
title('Respuesta del Empleo Informal a Shock de Productividad')
grid on

subplot(3,1,2)
plot(t, Lf_response, 'g-', 'LineWidth', 2)
hold on
plot([shock_period shock_period], [min(Lf_response) max(Lf_response)], 'r--')
xlabel('Período')
ylabel('L^f')
title('Respuesta del Empleo Formal')
grid on

subplot(3,1,3)
plot(t, r_response, 'm-', 'LineWidth', 2)
hold on
plot([shock_period shock_period], [min(r_response) max(r_response)], 'r--')
xlabel('Período')
ylabel('r')
title('Respuesta de la Tasa de Interés')
grid on

%% =========================================================================
% VALIDACIÓN CON DATOS EMPÍRICOS
% =========================================================================

fprintf('\n\nVALIDACIÓN CON DATOS EMPÍRICOS\n');
fprintf('================================\n\n');

% Datos de referencia (Schneider 2007, La Porta & Shleifer 2014)
empirical_data = struct();

countries = {'USA', 'España', 'México', 'Brasil', 'Colombia', 'Perú'};
informal_share = [0.08, 0.22, 0.30, 0.40, 0.38, 0.60];
gini_wealth = [0.80, 0.58, 0.77, 0.78, 0.77, 0.74];
tax_rate = [0.28, 0.35, 0.20, 0.34, 0.18, 0.17];

figure('Name', 'Validación Empírica', 'Position', [100 100 800 400])

subplot(1,2,1)
scatter(tax_rate, informal_share, 100, 'filled')
hold on
% Línea de tendencia
p = polyfit(tax_rate, informal_share, 1);
x_fit = linspace(min(tax_rate), max(tax_rate), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'r--', 'LineWidth', 2)
xlabel('Tasa Impositiva')
ylabel('Share Sector Informal')
title('Informalidad vs Impuestos')
grid on
% Etiquetas
for i = 1:length(countries)
    text(tax_rate(i)+0.01, informal_share(i)+0.02, countries{i}, 'FontSize', 8)
end

subplot(1,2,2)
scatter(informal_share, gini_wealth, 100, 'filled')
xlabel('Share Sector Informal')
ylabel('Gini de Riqueza')
title('Desigualdad vs Informalidad')
grid on
% Etiquetas
for i = 1:length(countries)
    text(informal_share(i)+0.01, gini_wealth(i)+0.01, countries{i}, 'FontSize', 8)
end

%% =========================================================================
% RECOMENDACIONES FINALES
% =========================================================================

fprintf('\n\n========== RECOMENDACIONES PARA IMPLEMENTACIÓN ==========\n\n');

fprintf('1. CALIBRACIÓN:\n');
fprintf('   - Usar datos de encuestas de hogares para calibrar shares sectoriales\n');
fprintf('   - La elasticidad de Frisch puede variar entre 0.5-2.0\n');
fprintf('   - Los parámetros ψᶠ y ψⁱ son cruciales para el tamaño sectorial\n\n');

fprintf('2. ASPECTOS COMPUTACIONALES:\n');
fprintf('   - El loop de punto fijo para wⁱ puede necesitar relajación\n');
fprintf('   - Considerar usar interpolación para acelerar convergencia\n');
fprintf('   - La grilla de activos puede necesitar más puntos cerca de ā\n\n');

fprintf('3. EXTENSIONES SUGERIDAS:\n');
fprintf('   - Añadir shocks de productividad sectorial\n');
fprintf('   - Incluir probabilidad de detección endógena\n');
fprintf('   - Considerar capital de trabajo en sector formal\n');
fprintf('   - Añadir costos de entrada/salida entre sectores\n\n');

fprintf('4. VALIDACIÓN:\n');
fprintf('   - Comparar momentos con datos de panel de empleo\n');
fprintf('   - Verificar transiciones sectoriales\n');
fprintf('   - Analizar correlación riqueza-formalidad\n\n');
