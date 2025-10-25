% AIYAGARI_DOS_SECTORES_ASSET_SUPPLY.m
% Versión ASSET SUPPLY del modelo con sectores formal e informal
% Solo resuelve el lado de los hogares para r exógena

clear all; clc; close all;

tic;

%% =========================================================================
% 1. PARÁMETROS DEL MODELO
% =========================================================================

% --- Parámetros de los Hogares ---
ga = 2;            % Coeficiente de aversión relativa al riesgo
rho = 0.05;        % Tasa de descuento subjetiva
frisch = 0.5;      % Elasticidad de Frisch

% --- Parámetros de Desutilidad del Trabajo ---
psi_f = 1.0;       % Desutilidad trabajo formal
psi_i = 0.9;       % Desutilidad trabajo informal

% --- Parámetros Fiscales ---
tau = 0.25;        % Tasa impositiva sobre trabajo formal

% --- Parámetros de Productividad EXÓGENOS (no hay firmas) ---
% En asset supply, tomamos salarios como dados
wf = 1.0;          % Salario formal EXÓGENO
wi = 0.7;          % Salario informal EXÓGENO

% --- Proceso Estocástico (Productividad Individual) ---
z1 = 0.5;          
z2 = 1.5;          
z = [z1,z2];       
la1 = 1/3;         
la2 = 1/3;         
la = [la1,la2];

%% =========================================================================
% 2. GRILLA DE TASAS DE INTERÉS A EVALUAR
% =========================================================================

% Vector de tasas de interés para calcular KS(r)
n_r = 10;
r_vec = linspace(0.01, 0.048, n_r);

% Para almacenar resultados
KS_vec = zeros(n_r, 1);
Lf_vec = zeros(n_r, 1);
Li_vec = zeros(n_r, 1);
gini_vec = zeros(n_r, 1);

%% =========================================================================
% 3. GRILLA DE ACTIVOS
% =========================================================================

I = 300;           % Menos puntos para velocidad
amin = 0;          
amax = 15;         

a = linspace(amin,amax,I)';
da = (amax-amin)/(I-1);

aa = [a,a];
zz = ones(I,1)*z;

%% =========================================================================
% 4. PARÁMETROS NUMÉRICOS
% =========================================================================

maxit = 100;
crit = 10^(-6);
Delta = 1000;

dVf = zeros(I,2);
dVb = zeros(I,2);

Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

%% =========================================================================
% 5. LOOP SOBRE TASAS DE INTERÉS (ASSET SUPPLY)
% =========================================================================

fprintf('CALCULANDO ASSET SUPPLY FUNCTION\n');
fprintf('=================================\n\n');

for i_r = 1:n_r
    
    r = r_vec(i_r);
    fprintf('Evaluando r = %.4f (%d de %d)\n', r, i_r, n_r);
    
    % --- 5.1. Calcular Oferta Laboral de Ahorro Cero ---
    lf0 = zeros(I,2);
    li0 = zeros(I,2);
    
    for j = 1:2
        for i = 1:I
            % Aproximación simple
            income_base = r * a(i);
            if income_base <= 0
                income_base = 0.01;
            end
            
            lf_temp = min(0.5, max(0, ((wf*z(j)*(1-tau))/psi_f)^frisch));
            li_temp = min(0.5, max(0, ((wi*z(j))/psi_i)^frisch));
            
            if lf_temp + li_temp > 1
                total = lf_temp + li_temp;
                lf_temp = lf_temp / total;
                li_temp = li_temp / total;
            end
            
            lf0(i,j) = lf_temp;
            li0(i,j) = li_temp;
        end
    end
    
    % --- 5.2. Conjetura Inicial para V ---
    if i_r == 1
        for j = 1:2
            income = wf*z(j)*lf0(:,j)*(1-tau) + wi*z(j)*li0(:,j) + r.*a;
            v0(:,j) = (income.^(1-ga))/(1-ga)/rho;
        end
    else
        v0 = V; % Warm start desde r anterior
    end
    
    v = v0;
    
    % --- 5.3. RESOLVER HJB ---
    for n = 1:maxit
        V = v;
        
        % Diferencias Finitas
        dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
        for j = 1:2
            income_max = wf*z(j)*(1-tau) + wi*z(j) + r*amax;
            dVf(I,j) = (income_max)^(-ga);
        end
        
        dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
        for j = 1:2
            income_min = wf*z(j)*lf0(1,j)*(1-tau) + wi*z(j)*li0(1,j) + r*amin;
            dVb(1,j) = (income_min)^(-ga);
        end
        
        % Esquema Upwind
        % Forward
        cf = dVf.^(-1/ga);
        lf_f = ((dVf.*wf.*zz*(1-tau))/psi_f).^frisch;
        li_f = ((dVf.*wi.*zz)/psi_i).^frisch;
        
        ltotal_f = lf_f + li_f;
        idx_exceed = ltotal_f > 1;
        lf_f(idx_exceed) = lf_f(idx_exceed) ./ ltotal_f(idx_exceed);
        li_f(idx_exceed) = li_f(idx_exceed) ./ ltotal_f(idx_exceed);
        
        ssf = wf*zz.*lf_f*(1-tau) + wi*zz.*li_f + r.*aa - cf;
        
        % Backward
        cb = dVb.^(-1/ga);
        lf_b = ((dVb.*wf.*zz*(1-tau))/psi_f).^frisch;
        li_b = ((dVb.*wi.*zz)/psi_i).^frisch;
        
        ltotal_b = lf_b + li_b;
        idx_exceed = ltotal_b > 1;
        lf_b(idx_exceed) = lf_b(idx_exceed) ./ ltotal_b(idx_exceed);
        li_b(idx_exceed) = li_b(idx_exceed) ./ ltotal_b(idx_exceed);
        
        ssb = wf*zz.*lf_b*(1-tau) + wi*zz.*li_b + r.*aa - cb;
        
        % Steady state
        c0 = wf*zz.*lf0*(1-tau) + wi*zz.*li0 + r.*aa;
        
        % Indicadores
        If = ssf > 0;
        Ib = ssb < 0;
        I0 = (1-If-Ib);
        
        % Políticas
        c = cf.*If + cb.*Ib + c0.*I0;
        lf = lf_f.*If + lf_b.*Ib + lf0.*I0;
        li = li_f.*If + li_b.*Ib + li0.*I0;
        
        % Utilidad
        u = c.^(1-ga)/(1-ga) - psi_f*lf.^(1+1/frisch)/(1+1/frisch) ...
            - psi_i*li.^(1+1/frisch)/(1+1/frisch);
        
        % Matriz A
        X = -min(ssb,0)/da;
        Y = -max(ssf,0)/da + min(ssb,0)/da;
        Z = max(ssf,0)/da;
        
        A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
        A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);
        
        A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
        
        % Sistema Implícito
        B = (1/Delta + rho)*speye(2*I) - A;
        
        u_stacked = [u(:,1);u(:,2)];
        V_stacked = [V(:,1);V(:,2)];
        
        b = u_stacked + V_stacked/Delta;
        V_stacked = B\b;
        
        V = [V_stacked(1:I),V_stacked(I+1:2*I)];
        
        % Convergencia
        Vchange = V - v;
        v = V;
        
        dist(n) = max(max(abs(Vchange)));
        if dist(n) < crit
            break
        end
    end
    
    % --- 5.4. FOKKER-PLANCK ---
    AT = A';
    b = zeros(2*I,1);
    
    i_fix = 1;
    b(i_fix) = .1;
    row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
    AT(i_fix,:) = row;
    
    gg = AT\b;
    g_sum = gg'*ones(2*I,1)*da;
    gg = gg./g_sum;
    
    g = [gg(1:I),gg(I+1:2*I)];
    
    % --- 5.5. CALCULAR Y ALMACENAR AGREGADOS ---
    
    % Capital Supply
    KS_vec(i_r) = g(:,1)'*a*da + g(:,2)'*a*da;
    
    % Trabajo agregado por sector
    Lf_vec(i_r) = 0;
    Li_vec(i_r) = 0;
    for j = 1:2
        Lf_vec(i_r) = Lf_vec(i_r) + g(:,j)'*lf(:,j)*da;
        Li_vec(i_r) = Li_vec(i_r) + g(:,j)'*li(:,j)*da;
    end
    
    % Gini
    g_total = g(:,1) + g(:,2);
    [a_sorted, idx] = sort(a);
    g_sorted = g_total(idx);
    g_cum = cumsum(g_sorted)*da;
    a_cum = cumsum(a_sorted.*g_sorted)*da;
    a_total = a_cum(end);
    gini_vec(i_r) = 1 - 2*sum(a_cum)*da/a_total;
    
    fprintf('  KS = %.4f, Lf = %.4f, Li = %.4f, Gini = %.4f\n', ...
        KS_vec(i_r), Lf_vec(i_r), Li_vec(i_r), gini_vec(i_r));
    
    % Guardar para última r (para gráficos)
    if i_r == n_r
        g_final = g;
        lf_final = lf;
        li_final = li;
        c_final = c;
        s_final = wf*zz.*lf*(1-tau) + wi*zz.*li + r.*aa - c;
    end
end

toc;

%% =========================================================================
% 6. GRÁFICOS DE RESULTADOS
% =========================================================================

% Figura 1: Asset Supply Function
figure(1)
subplot(2,2,1)
plot(r_vec*100, KS_vec, 'b-o', 'LineWidth', 2)
xlabel('Tasa de Interés r (%)')
ylabel('Capital Supply KS')
title('Asset Supply Function')
grid on

subplot(2,2,2)
plot(r_vec*100, Lf_vec, 'r-o', 'LineWidth', 2)
hold on
plot(r_vec*100, Li_vec, 'g-o', 'LineWidth', 2)
xlabel('Tasa de Interés r (%)')
ylabel('Empleo')
title('Empleo por Sector')
legend('Formal', 'Informal', 'Location', 'best')
grid on

subplot(2,2,3)
plot(r_vec*100, Li_vec./(Lf_vec+Li_vec)*100, 'm-o', 'LineWidth', 2)
xlabel('Tasa de Interés r (%)')
ylabel('Share Informal (%)')
title('Participación del Sector Informal')
grid on

subplot(2,2,4)
plot(r_vec*100, gini_vec, 'k-o', 'LineWidth', 2)
xlabel('Tasa de Interés r (%)')
ylabel('Índice de Gini')
title('Desigualdad')
grid on

% Figura 2: Políticas para última r
figure(2)
amax_plot = 5;

subplot(2,2,1)
plot(a, s_final, 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Ahorro s(a,z)')
title(sprintf('Política de Ahorro (r = %.3f)', r))
legend('z_1', 'z_2', 'Location', 'best')
xlim([0 amax_plot])
grid on

subplot(2,2,2)
plot(a, c_final, 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Consumo c(a,z)')
title('Política de Consumo')
legend('z_1', 'z_2', 'Location', 'best')
xlim([0 amax_plot])
grid on

subplot(2,2,3)
plot(a, lf_final, 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Horas Formales')
title('Trabajo Formal')
legend('z_1', 'z_2', 'Location', 'best')
xlim([0 amax_plot])
grid on

subplot(2,2,4)
plot(a, li_final, 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Horas Informales')
title('Trabajo Informal')
legend('z_1', 'z_2', 'Location', 'best')
xlim([0 amax_plot])
grid on

% Figura 3: Distribución
figure(3)
subplot(1,2,1)
plot(a, g_final, 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Densidad g(a,z)')
title(sprintf('Distribución Estacionaria (r = %.3f)', r))
legend('z_1', 'z_2', 'Location', 'best')
xlim([0 amax_plot])
grid on

subplot(1,2,2)
bar_data = [sum(g_final(:,1))*da, sum(g_final(:,2))*da];
bar(1:2, bar_data)
set(gca, 'XTickLabel', {'z_1 (baja)', 'z_2 (alta)'})
ylabel('Masa de Probabilidad')
title('Distribución de Productividad')
grid on

%% =========================================================================
% 7. EXPORTAR RESULTADOS
% =========================================================================

% Guardar asset supply function
results.r_vec = r_vec;
results.KS_vec = KS_vec;
results.Lf_vec = Lf_vec;
results.Li_vec = Li_vec;
results.gini_vec = gini_vec;
results.share_informal = Li_vec./(Lf_vec+Li_vec);

% Guardar en archivo
save('asset_supply_results.mat', 'results', 'g_final', 'lf_final', 'li_final');

fprintf('\n========== RESUMEN ASSET SUPPLY ==========\n');
fprintf('Rango de r evaluado: [%.3f, %.3f]\n', min(r_vec), max(r_vec));
fprintf('Rango de KS: [%.2f, %.2f]\n', min(KS_vec), max(KS_vec));
fprintf('Elasticidad KS a r (promedio): %.3f\n', ...
    mean(diff(log(KS_vec))./diff(log(r_vec))));
fprintf('Share informal promedio: %.2f%%\n', mean(results.share_informal)*100);
fprintf('\nResultados guardados en: asset_supply_results.mat\n');

%% =========================================================================
% 8. ENCONTRAR EQUILIBRIO (OPCIONAL)
% =========================================================================

% Si tuviéramos una función de demanda de capital KD(r), 
% podríamos encontrar el equilibrio

% Ejemplo: supongamos KD = A*(r+d)^(-epsilon)
A_demand = 10;
d = 0.05;
epsilon = 2;
KD_vec = A_demand * (r_vec + d).^(-epsilon);

% Encontrar intersección
figure(4)
plot(r_vec*100, KS_vec, 'b-o', 'LineWidth', 2)
hold on
plot(r_vec*100, KD_vec, 'r--', 'LineWidth', 2)
xlabel('Tasa de Interés r (%)')
ylabel('Capital')
title('Equilibrio en Mercado de Capital')
legend('Supply KS(r)', 'Demand KD(r)', 'Location', 'best')
grid on

% Encontrar r de equilibrio (donde KS = KD)
[~, idx_eq] = min(abs(KS_vec - KD_vec));
r_eq = r_vec(idx_eq);
K_eq = KS_vec(idx_eq);

plot(r_eq*100, K_eq, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
text(r_eq*100+0.2, K_eq, sprintf('Equilibrio\nr* = %.3f\nK* = %.2f', r_eq, K_eq))

fprintf('\n========== EQUILIBRIO (SI KD FUERA CONOCIDA) ==========\n');
fprintf('r* = %.4f\n', r_eq);
fprintf('K* = %.4f\n', K_eq);
fprintf('Share informal en equilibrio: %.2f%%\n', ...
    results.share_informal(idx_eq)*100);
