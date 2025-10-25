% AIYAGARI_DOS_SECTORES.m
% Modelo de Aiyagari con trabajo endógeno en sectores formal e informal
% Basado en Moll et al. + Restrepo-Echavarria

clear all; clc; close all;

tic; % Inicia el cronómetro

%% =========================================================================
% 1. PARÁMETROS DEL MODELO
% =========================================================================

% --- Parámetros de los Hogares ---
ga = 2;            % Coeficiente de aversión relativa al riesgo (γ)
rho = 0.05;        % Tasa de descuento subjetiva
d = 0.05;          % Tasa de depreciación del capital (δ)
frisch = 0.5;      % Elasticidad de Frisch (φ)

% --- Parámetros de Desutilidad del Trabajo ---
psi_f = 1.0;       % Desutilidad trabajo formal (ψᶠ)
psi_i = 0.9;       % Desutilidad trabajo informal (ψⁱ) - Menor = informal más atractivo

% --- Parámetros Fiscales ---
tau = 0.25;        % Tasa impositiva sobre trabajo formal

% --- Parámetros de las Firmas ---
% Sector Formal (Cobb-Douglas con capital)
al_f = 1/3;        % Participación del capital en formal (αᶠ)
Af = 0.15;         % Productividad total formal

% Sector Informal (Solo trabajo, retornos decrecientes)
al_i = 0.75;       % Exponente trabajo informal (αⁱ < 1)
Ai = 0.08;         % Productividad informal

% --- Parámetros del Proceso Estocástico (Productividad Individual) ---
z1 = 0.5;          % Productividad baja
z2 = 1.5;          % Productividad alta
z = [z1,z2];       % Vector de estados
la1 = 1/3;         % Tasa de Poisson z1 -> z2
la2 = 1/3;         % Tasa de Poisson z2 -> z1
la = [la1,la2];

% Productividad promedio
z_ave = (z1*la2 + z2*la1)/(la1 + la2);

%% =========================================================================
% 2. GRILLA DE ACTIVOS
% =========================================================================

I = 500;           % Número de puntos en la grilla (reducido para velocidad)
amin = 0;          % Límite inferior (restricción de no-endeudamiento)
amax = 15;         % Límite superior de activos

a = linspace(amin,amax,I)'; % Vector columna de la grilla
da = (amax-amin)/(I-1);     % Distancia entre puntos

% Matrices auxiliares
aa = [a,a];        % Repite la grilla para ambos estados de z
zz = ones(I,1)*z;  % Matriz (I x 2) con valores de z

%% =========================================================================
% 3. PARÁMETROS DE ITERACIÓN NUMÉRICA
% =========================================================================

% --- Parámetros HJB ---
maxit = 100;       % Máximo de iteraciones para HJB
crit = 10^(-6);    % Criterio de convergencia HJB
Delta = 1000;      % Tamaño del paso para método implícito

% Inicialización de matrices
dVf = zeros(I,2);  % Derivada forward
dVb = zeros(I,2);  % Derivada backward
c = zeros(I,2);    % Política de consumo
lf = zeros(I,2);   % Política trabajo formal
li = zeros(I,2);   % Política trabajo informal

% Matriz de saltos entre estados
Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

% --- Parámetros Equilibrio General ---
Ir = 30;           % Máximas iteraciones para encontrar r
crit_S = 10^(-5);  % Criterio convergencia mercado capital
crit_L = 10^(-4);  % Criterio convergencia mercado trabajo informal

% --- Parámetros de Bisección ---
rmin = 0.01;       % Límite inferior para r
rmax = 0.99*rho;   % Límite superior para r
r = 0.04;          % Conjetura inicial

%% =========================================================================
% 4. BUCLE EXTERNO: EQUILIBRIO GENERAL
% =========================================================================

for ir = 1:Ir
    
    % Almacena historial
    r_r(ir) = r;
    rmin_r(ir) = rmin;
    rmax_r(ir) = rmax;
    
    fprintf('\n========== Iteración Equilibrio %d ==========\n', ir);
    fprintf('Tasa de interés r = %.4f\n', r);
    
    % --- 4.1. Resolver Precios del Sector Formal (dado r) ---
    % FOC Capital Formal: r + d = al_f * Af * (Lf/Kf)^(1-al_f)
    % Por ahora asumimos Lf = z_ave para inicializar
    Lf_guess = z_ave * 0.7; % 70% en formal inicialmente
    KD = (al_f * Af / (r + d))^(1/(1-al_f)) * Lf_guess;
    wf = (1-al_f) * Af * (KD/Lf_guess)^al_f; % Salario formal
    
    % --- 4.2. Loop de Punto Fijo para Sector Informal ---
    % El salario informal depende del trabajo informal agregado
    
    Li_old = z_ave * 0.3; % 30% en informal inicialmente
    for iter_wi = 1:20
        % Salario informal desde FOC: wi = al_i * Ai * Li^(al_i-1)
        wi = al_i * Ai * Li_old^(al_i-1);
        
        fprintf('  Iteración salario informal %d: wi = %.4f\n', iter_wi, wi);
        
        % --- 4.3. Calcular Oferta Laboral de Ahorro Cero ---
        % Necesitamos resolver para (lf0, li0) cuando s=0
        fprintf('  Calculando oferta laboral de ahorro cero...\n');
        
        lf0 = zeros(I,2);
        li0 = zeros(I,2);
        
        for j = 1:2  % Para cada estado de productividad
            for i = 1:I  % Para cada nivel de activos
                % Resolver sistema: s = 0 y FOCs de trabajo
                % wf*z*lf*(1-tau) + wi*z*li + r*a = c
                % psi_f * lf^(1/frisch) = c^(-ga) * wf*z*(1-tau)
                % psi_i * li^(1/frisch) = c^(-ga) * wi*z
                
                % Simplificación: usar aproximación
                income_base = r * a(i);
                if income_base <= 0
                    income_base = 0.01; % Evitar problemas numéricos
                end
                
                % Aproximación inicial
                lf_temp = min(0.5, max(0, ((wf*z(j)*(1-tau))/psi_f)^frisch));
                li_temp = min(0.5, max(0, ((wi*z(j))/psi_i)^frisch));
                
                % Normalizar si suma > 1
                if lf_temp + li_temp > 1
                    total = lf_temp + li_temp;
                    lf_temp = lf_temp / total;
                    li_temp = li_temp / total;
                end
                
                lf0(i,j) = lf_temp;
                li0(i,j) = li_temp;
            end
        end
        
        % --- 4.4. Conjetura Inicial para Función de Valor ---
        if ir == 1 && iter_wi == 1
            % Primera iteración: usar aproximación
            for j = 1:2
                income = wf*z(j)*lf0(:,j)*(1-tau) + wi*z(j)*li0(:,j) + r.*a;
                v0(:,j) = (income.^(1-ga))/(1-ga)/rho;
            end
        else
            % Usar valor previo (warm start)
            if iter_wi == 1
                v0 = V_r(:,:,ir-1);
            else
                v0 = V;
            end
        end
        
        v = v0;
        
        % --- 4.5. BUCLE INTERNO: RESOLVER HJB ---
        for n = 1:maxit
            V = v;
            
            % Derivadas (Diferencias Finitas)
            % Forward
            dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
            for j = 1:2
                income_max = wf*z(j)*(1-tau) + wi*z(j) + r*amax;
                dVf(I,j) = (income_max)^(-ga);
            end
            
            % Backward
            dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
            for j = 1:2
                income_min = wf*z(j)*lf0(1,j)*(1-tau) + wi*z(j)*li0(1,j) + r*amin;
                dVb(1,j) = (income_min)^(-ga);
            end
            
            % --- ESQUEMA UPWIND ---
            % Forward
            cf = dVf.^(-1/ga);
            lf_f = ((dVf.*wf.*zz*(1-tau))/psi_f).^frisch;
            li_f = ((dVf.*wi.*zz)/psi_i).^frisch;
            
            % Restringir horas totales
            ltotal_f = lf_f + li_f;
            idx_exceed = ltotal_f > 1;
            lf_f(idx_exceed) = lf_f(idx_exceed) ./ ltotal_f(idx_exceed);
            li_f(idx_exceed) = li_f(idx_exceed) ./ ltotal_f(idx_exceed);
            
            ssf = wf*zz.*lf_f*(1-tau) + wi*zz.*li_f + r.*aa - cf;
            
            % Backward
            cb = dVb.^(-1/ga);
            lf_b = ((dVb.*wf.*zz*(1-tau))/psi_f).^frisch;
            li_b = ((dVb.*wi.*zz)/psi_i).^frisch;
            
            % Restringir horas totales
            ltotal_b = lf_b + li_b;
            idx_exceed = ltotal_b > 1;
            lf_b(idx_exceed) = lf_b(idx_exceed) ./ ltotal_b(idx_exceed);
            li_b(idx_exceed) = li_b(idx_exceed) ./ ltotal_b(idx_exceed);
            
            ssb = wf*zz.*lf_b*(1-tau) + wi*zz.*li_b + r.*aa - cb;
            
            % En steady state (s=0)
            c0 = wf*zz.*lf0*(1-tau) + wi*zz.*li0 + r.*aa;
            
            % Indicadores upwind
            If = ssf > 0;
            Ib = ssb < 0;
            I0 = (1-If-Ib);
            
            % Políticas finales
            c = cf.*If + cb.*Ib + c0.*I0;
            lf = lf_f.*If + lf_b.*Ib + lf0.*I0;
            li = li_f.*If + li_b.*Ib + li0.*I0;
            
            % Utilidad
            u = c.^(1-ga)/(1-ga) - psi_f*lf.^(1+1/frisch)/(1+1/frisch) ...
                - psi_i*li.^(1+1/frisch)/(1+1/frisch);
            
            % --- MATRIZ DE TRANSICIÓN ---
            X = -min(ssb,0)/da;
            Y = -max(ssf,0)/da + min(ssb,0)/da;
            Z = max(ssf,0)/da;
            
            A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);
            
            A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
            
            % --- SISTEMA IMPLÍCITO ---
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
                fprintf('    HJB convergió en %d iteraciones\n', n);
                break
            end
        end
        
        % --- 4.6. FOKKER-PLANCK ---
        AT = A';
        b = zeros(2*I,1);
        
        % Fix: primera condición
        i_fix = 1;
        b(i_fix) = .1;
        row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
        AT(i_fix,:) = row;
        
        gg = AT\b;
        g_sum = gg'*ones(2*I,1)*da;
        gg = gg./g_sum;
        
        g = [gg(1:I),gg(I+1:2*I)];
        
        % --- 4.7. Calcular Agregados ---
        % Trabajo formal agregado
        Lf_new = 0;
        for j = 1:2
            Lf_new = Lf_new + g(:,j)'*lf(:,j)*da;
        end
        
        % Trabajo informal agregado
        Li_new = 0;
        for j = 1:2
            Li_new = Li_new + g(:,j)'*li(:,j)*da;
        end
        
        fprintf('    Lf = %.4f, Li = %.4f\n', Lf_new, Li_new);
        
        % Verificar convergencia del trabajo informal
        if abs(Li_new - Li_old) < crit_L
            fprintf('  Mercado trabajo informal convergió\n');
            Lf = Lf_new;
            Li = Li_new;
            break
        end
        
        % Actualizar
        Li_old = 0.5*Li_old + 0.5*Li_new; % Promedio para estabilidad
        
        % Actualizar capital y salario formal
        KD = (al_f * Af / (r + d))^(1/(1-al_f)) * Lf_new;
        wf = (1-al_f) * Af * (KD/Lf_new)^al_f;
    end
    
    % Almacenar resultados de esta iteración
    g_r(:,:,ir) = g;
    V_r(:,:,ir) = V;
    lf_r(:,:,ir) = lf;
    li_r(:,:,ir) = li;
    
    % --- 4.8. EQUILIBRIO DE CAPITAL ---
    % Capital ofrecido por hogares
    KS = g(:,1)'*a*da + g(:,2)'*a*da;
    
    % Exceso de oferta
    S(ir) = KS - KD;
    
    fprintf('KS = %.4f, KD = %.4f, Exceso = %.4f\n', KS, KD, S(ir));
    
    % Bisección
    if S(ir) > crit_S
        fprintf('Exceso de Oferta de Capital\n');
        rmax = r;
        r = 0.5*(r+rmin);
    elseif S(ir) < -crit_S
        fprintf('Exceso de Demanda de Capital\n');
        rmin = r;
        r = 0.5*(r+rmax);
    else
        fprintf('\n¡EQUILIBRIO ENCONTRADO!\n');
        fprintf('r* = %.4f, Lf = %.4f, Li = %.4f\n', r, Lf, Li);
        break
    end
end

toc;

%% =========================================================================
% 5. GRÁFICOS DE RESULTADOS
% =========================================================================

amax_plot = 5;

% Figura 1: Políticas de Ahorro
figure(1)
subplot(1,2,1)
plot(a, wf*zz.*lf_r(:,:,ir)*(1-tau) + wi*zz.*li_r(:,:,ir) + r.*aa - c, 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Ahorro s(a,z)')
title('Política de Ahorro')
legend('z_1 (baja)', 'z_2 (alta)', 'Location', 'best')
xlim([0 amax_plot])
grid on

subplot(1,2,2)
plot(a, g_r(:,:,ir), 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Densidad g(a,z)')
title('Distribución Estacionaria')
legend('z_1 (baja)', 'z_2 (alta)', 'Location', 'best')
xlim([0 amax_plot])
grid on

% Figura 2: Asignación Laboral
figure(2)
subplot(2,2,1)
plot(a, lf_r(:,:,ir), 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Horas Formales l^f(a,z)')
title('Trabajo Formal')
legend('z_1 (baja)', 'z_2 (alta)', 'Location', 'best')
xlim([0 amax_plot])
grid on

subplot(2,2,2)
plot(a, li_r(:,:,ir), 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Horas Informales l^i(a,z)')
title('Trabajo Informal')
legend('z_1 (baja)', 'z_2 (alta)', 'Location', 'best')
xlim([0 amax_plot])
grid on

subplot(2,2,3)
plot(a, lf_r(:,:,ir) + li_r(:,:,ir), 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Horas Totales l^f + l^i')
title('Trabajo Total')
legend('z_1 (baja)', 'z_2 (alta)', 'Location', 'best')
xlim([0 amax_plot])
grid on

subplot(2,2,4)
plot(a, lf_r(:,:,ir)./(lf_r(:,:,ir) + li_r(:,:,ir)), 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Share Formal l^f/(l^f+l^i)')
title('Proporción en Sector Formal')
legend('z_1 (baja)', 'z_2 (alta)', 'Location', 'best')
xlim([0 amax_plot])
ylim([0 1])
grid on

% Figura 3: Consumo
figure(3)
plot(a, c, 'LineWidth', 2)
xlabel('Riqueza (a)')
ylabel('Consumo c(a,z)')
title('Política de Consumo')
legend('z_1 (baja)', 'z_2 (alta)', 'Location', 'best')
xlim([0 amax_plot])
grid on

% Figura 4: Estadísticas Agregadas
figure(4)
% Calcular share informal por deciles de riqueza
n_deciles = 10;
for d = 1:n_deciles
    a_min_d = (d-1)*amax/n_deciles;
    a_max_d = d*amax/n_deciles;
    idx_d = (a >= a_min_d) & (a < a_max_d);
    
    % Peso de cada decil
    mass_d = sum(sum(g(idx_d,:)))*da;
    
    if mass_d > 0
        % Horas promedio en informal para este decil
        li_avg = 0;
        lf_avg = 0;
        for j = 1:2
            li_avg = li_avg + sum(g(idx_d,j).*li_r(idx_d,j,ir))*da/mass_d;
            lf_avg = lf_avg + sum(g(idx_d,j).*lf_r(idx_d,j,ir))*da/mass_d;
        end
        share_informal(d) = li_avg/(li_avg + lf_avg);
    else
        share_informal(d) = 0;
    end
end

bar(1:n_deciles, share_informal)
xlabel('Decil de Riqueza')
ylabel('Share de Trabajo Informal')
title('Informalidad por Decil de Riqueza')
grid on

%% =========================================================================
% 6. ESTADÍSTICAS RESUMEN
% =========================================================================

fprintf('\n========== ESTADÍSTICAS DEL EQUILIBRIO ==========\n');
fprintf('Tasa de interés de equilibrio: r* = %.4f\n', r);
fprintf('Salario formal: wf = %.4f\n', wf);
fprintf('Salario informal: wi = %.4f\n', wi);
fprintf('Capital: K = %.4f\n', KS);
fprintf('Trabajo formal agregado: Lf = %.4f\n', Lf);
fprintf('Trabajo informal agregado: Li = %.4f\n', Li);
fprintf('Share informal en empleo: %.2f%%\n', 100*Li/(Lf+Li));

% Calcular producto
Yf = Af * KD^al_f * Lf^(1-al_f);
Yi = Ai * Li^al_i;
Y_total = Yf + Yi;

fprintf('Producto formal: Yf = %.4f\n', Yf);
fprintf('Producto informal: Yi = %.4f\n', Yi);
fprintf('Share informal en producto: %.2f%%\n', 100*Yi/Y_total);

% Gini de riqueza
[a_sorted, idx] = sort(a);
g_total = g(:,1) + g(:,2);
g_sorted = g_total(idx);
g_cum = cumsum(g_sorted)*da;
a_cum = cumsum(a_sorted.*g_sorted)*da;
a_total = a_cum(end);
gini = 1 - 2*sum(a_cum)*da/a_total;
fprintf('Índice de Gini (riqueza): %.4f\n', gini);
