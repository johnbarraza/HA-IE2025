function [lf, li, c] = lab_solve_dos_sectores(a, z, wf, wi, r, ga, frisch, psi_f, psi_i, tau)
% LAB_SOLVE_DOS_SECTORES - Resuelve oferta laboral óptima con ahorro cero (s=0)
% para modelo con sectores formal e informal. Maneja restricciones.

% Opciones para el solver
options_fsolve = optimset('Display', 'off', 'TolFun', 1e-9, 'TolX', 1e-9);
options_fzero = optimset('Display','off', 'TolFun', 1e-9, 'TolX', 1e-9); % Esta ya estaba bien

% --- Intento 1: Resolver el sistema de FOCs sin imponer l^f + l^i <= 1 ---
fun = @(x) system_focs(x, a, z, wf, wi, r, ga, frisch, psi_f, psi_i, tau);

% Conjetura inicial (mejorada)
wf_eff = wf*z*(1-tau);
wi_eff = wi*z;
c_init_guess = max(r*a + 0.5*wf_eff + 0.5*wi_eff, 0.1); % Aprox c = ra + Y_medio
lf_init = max(0, min(1, (c_init_guess^(-ga) * wf_eff / psi_f)^frisch));
li_init = max(0, min(1, (c_init_guess^(-ga) * wi_eff / psi_i)^frisch));
x0 = [lf_init; li_init; c_init_guess];

[x_sol, fval, exitflag] = fsolve(fun, x0, options_fsolve);

lf_sol = x_sol(1);
li_sol = x_sol(2);
c_sol = x_sol(3);

% --- Comprobar validez de la solución de fsolve ---
valid_fsolve = (exitflag > 0) && (c_sol > 1e-9) && (lf_sol >= -1e-9) && (li_sol >= -1e-9) && ~any(isnan(x_sol));

if valid_fsolve && (lf_sol + li_sol <= 1 + 1e-9) % Solución interior o justo en la frontera
    lf = max(0, lf_sol);
    li = max(0, li_sol);
    % Asegurar que la suma no exceda 1 por errores numéricos
    if lf + li > 1
        total = lf + li;
        lf = lf / total;
        li = li / total;
    end
    c = max(1e-9, c_sol);
    % Recalcular c por consistencia si hubo ajuste
    c = wf*z*lf*(1-tau) + wi*z*li + r*a;
    
else % fsolve falló o la solución está fuera de la frontera l^f + l^i = 1
    
    % --- Intento 2: Resolver sobre la frontera l^f + l^i = 1 ---
    % FOC_f / FOC_i => (lf / (1-lf))^(1/frisch) = (psi_i * wf * (1-tau)) / (psi_f * wi)

   lf_bound = NaN; % Inicializa
    util_bound = -Inf; % Inicializa

    % Calcular ratio_target de forma segura
    denominator_ratio = psi_f * wi * z_val;
    if abs(denominator_ratio) < 1e-20
        if abs(psi_i * wf * z_val * (1-tau)) < 1e-20
             ratio_target = NaN; % Caso 0/0
        else
             ratio_target = Inf; % Caso num/0
        end
    else
        ratio_target = (psi_i * wf * z_val * (1-tau)) / denominator_ratio;
    end

    % --- Manejo de casos extremos de ratio_target ---
    if isnan(ratio_target) || ratio_target < 0
        warning('Ratio target inválido (NaN o negativo). Asignando fallback.');
        if wf*(1-tau) > wi && wf*(1-tau) > 0
            lf_bound = 1;
        elseif wi > wf*(1-tau) && wi > 0
            lf_bound = 0;
        else
             lf_bound = 0;
        end
        
    elseif isinf(ratio_target) || ratio_target > 1e30
        lf_bound = 1;
    elseif ratio_target < 1e-30
        lf_bound = 0;
    else
        % --- ratio_target es finito, positivo y razonable: intentar evaluar ---
        fun_boundary = @(lf_b) (max(lf_b, 1e-10) / max(1 - lf_b, 1e-10))^(1/frisch) - ratio_target;

        % --- Evaluar en los extremos con manejo de errores y magnitud ---
        val_low = NaN; val_high = NaN;
        MAX_ACCEPTABLE_VALUE = 1e17; % Límite antes de considerar infinito numérico

        try
            val_low_temp = fun_boundary(1e-9);
            % Chequea finitud Y magnitud razonable
            if isfinite(val_low_temp) && abs(val_low_temp) < MAX_ACCEPTABLE_VALUE
                val_low = val_low_temp;
            else
                warning('Valor extremo o no finito en extremo inferior de fun_boundary.');
            end
        catch ME_low
            warning('Error evaluando fun_boundary en extremo inferior: %s', ME_low.message);
        end

        try
            val_high_temp = fun_boundary(1-1e-9);
            % Chequea finitud Y magnitud razonable
            if isfinite(val_high_temp) && abs(val_high_temp) < MAX_ACCEPTABLE_VALUE
                val_high = val_high_temp;
            else
                 warning('Valor extremo o no finito en extremo superior de fun_boundary.');
                 % Si val_high explotó, lf_bound probablemente es 1
                 if isnan(lf_bound) % Solo asigna si no se determinó antes
                     lf_bound = 1;
                 end
            end
        catch ME_high
            warning('Error evaluando fun_boundary en extremo superior: %s', ME_high.message);
        end

        % --- Decidir basado en los valores (si lf_bound no se asignó ya) ---
        if isnan(lf_bound) % Si no fue asignado por ratio_target extremo o por val_high extremo

            if isnan(val_low) || isnan(val_high) % Si alguna evaluación falló
                 warning('Valores no válidos en extremos de fun_boundary. Asignando según ratio_target.');
                 if ratio_target > 1
                     lf_bound = 1;
                 else
                     lf_bound = 0;
                 end
            else
                % --- Los valores son finitos y razonables, proceder ---
                if val_low * val_high < 0 && abs(val_low) > 1e-12 && abs(val_high) > 1e-12
                    % Hay cambio de signo claro, intentar fzero
                    [lf_bound_sol, ~, flag_fzero_bnd] = fzero(fun_boundary, [1e-9, 1-1e-9], options_fzero);
                    if flag_fzero_bnd == 1
                        lf_bound = lf_bound_sol;
                    else
                         warning('fzero falló en frontera, usando aproximación de borde.');
                         if abs(val_low) < abs(val_high)
                             lf_bound = 1e-9;
                         else
                             lf_bound = 1-1e-9;
                         end
                    end
                elseif abs(val_low) <= 1e-12
                     lf_bound = 0;
                elseif abs(val_high) <= 1e-12
                     lf_bound = 1;
                elseif val_low > 0
                     lf_bound = 0;
                elseif val_high < 0
                     lf_bound = 1;
                else
                     warning('Caso residual fun_boundary.');
                     if abs(val_low) < abs(val_high)
                          lf_bound = 0;
                     else
                          lf_bound = 1;
                     end
                end
            end
        end % Fin del if isnan(lf_bound)
    end % Fin del else (ratio_target razonable)

    % --- Calcular li_bound, c_bound, util_bound ---
    lf_bound = max(0, min(1, lf_bound));
    li_bound = 1 - lf_bound;
    c_bound = wf*z_val*lf_bound*(1-tau) + wi*z_val*li_bound + r*a; % Usa z_val
    c_bound = max(1e-9, c_bound);

    util_bound = -Inf;
    if c_bound > 1e-9 && isfinite(c_bound)
        util_val_f = - psi_f*lf_bound^(1+1/frisch)/(1+1/frisch);
        util_val_i = - psi_i*li_bound^(1+1/frisch)/(1+1/frisch);
        if isnan(util_val_f), util_val_f = 0; end
        if isnan(util_val_i), util_val_i = 0; end
        util_bound = c_bound^(1-ga)/(1-ga) + util_val_f + util_val_i;
        if isnan(util_bound), util_bound = -Inf; end
    end
    
    % --- Intento 3: Comprobar Esquinas ---
    
    % Esquina 1: lf=0, li=?
    [li_c1, c_c1, util_c1, flag_c1] = solve_corner(a, z, wf, wi, r, ga, frisch, psi_f, psi_i, tau, 'informal');
    
    % Esquina 2: li=0, lf=?
    [lf_c2, c_c2, util_c2, flag_c2] = solve_corner(a, z, wf, wi, r, ga, frisch, psi_f, psi_i, tau, 'formal');
    
    % Esquina 0: lf=0, li=0
    c_c0 = r*a;
    util_c0 = -Inf; % Asigna utilidad muy baja si c<=0 o si no trabaja
    if c_c0 > 0
        util_c0 = c_c0^(1-ga)/(1-ga); % Desutilidad es 0
    end
    
    % --- Comparar Solución de Frontera y Esquinas ---
    best_util = -Inf;
    lf = 0; li = 0; c = max(r*a, 1e-9); % Default a esquina 0
    
    % Comprobar Frontera
    if c_bound > 0 && ~isnan(util_bound)
        best_util = util_bound;
        lf = lf_bound;
        li = li_bound;
        c = c_bound;
    end
    
    % Comprobar Esquina 1 (Informal solo)
    if flag_c1 && util_c1 > best_util
        best_util = util_c1;
        lf = 0;
        li = li_c1;
        c = c_c1;
    end
    
    % Comprobar Esquina 2 (Formal solo)
    if flag_c2 && util_c2 > best_util
        % best_util = util_c2; % No hace falta reasignar best_util
        lf = lf_c2;
        li = 0;
        c = c_c2;
    end

    % Comprobar Esquina 0 (No trabajar)
    if util_c0 > best_util
       % best_util = util_c0; % No hace falta reasignar best_util
        lf = 0;
        li = 0;
        c = max(r*a, 1e-9);
    end
    
end % Fin del manejo si fsolve falló o violó restricción

% Asegurar positividad final del consumo
c = max(c, 1e-9);

end

%% =========================================================================
% Función Auxiliar: Sistema de FOCs (para fsolve, caso no restringido)
% =========================================================================
function F = system_focs(x, a, z, wf, wi, r, ga, frisch, psi_f, psi_i, tau)
    lf = x(1);
    li = x(2);
    c = x(3);
    
    % Prevenir valores negativos/cero problemáticos durante la iteración
    lf_eval = max(lf, 1e-10);
    li_eval = max(li, 1e-10);
    c_eval = max(c, 1e-10);
    
    F = zeros(3,1);
    
    % Ecuación 1: Restricción presupuestaria con s=0
    F(1) = c_eval - (wf*z*lf*(1-tau) + wi*z*li + r*a); % Usa lf, li originales aquí
    
    % Ecuación 2: FOC trabajo formal (MRS = Salario Neto Marginal)
    F(2) = psi_f * lf_eval^(1/frisch) - c_eval^(-ga) * wf*z*(1-tau);
    
    % Ecuación 3: FOC trabajo informal (MRS = Salario Neto Marginal)
    F(3) = psi_i * li_eval^(1/frisch) - c_eval^(-ga) * wi*z;
    
    % Penalizar si horas negativas o consumo negativo (para guiar a fsolve)
    penalty = 0;
    if lf < 0, penalty = penalty + 100*lf^2; end
    if li < 0, penalty = penalty + 100*li^2; end
    if c < 0, penalty = penalty + 100*c^2; end
    F = F + penalty;
    
end

%% =========================================================================
% Función Auxiliar: Resolver Esquinas (lf=0 o li=0)
% =========================================================================
function [h_opt, c_opt, util_opt, flag] = solve_corner(a, z, wf, wi, r, ga, frisch, psi_f, psi_i, tau, sector)
    
    options_fzero = optimset('Display','off', 'TolFun', 1e-9, 'TolX', 1e-9);
    h_opt = 0; c_opt = 0; util_opt = -Inf; flag = false; % Inicializa

    if strcmp(sector, 'formal')
        psi = psi_f;
        wage_eff = wf*z*(1-tau);
        other_wage_eff = wi*z; % Salario del sector excluido
        other_psi = psi_i;
    else % informal
        psi = psi_i;
        wage_eff = wi*z;
        other_wage_eff = wf*z*(1-tau);
        other_psi = psi_f;
    end
    
    if wage_eff <= 0 % Si el salario efectivo es no positivo, la solución es 0
        h_opt = 0;
        c_opt = r*a; % Solo ingreso de capital
    else
        % Función para FOC de esquina: MRS = Salario efectivo
        % psi * h^(1/frisch) = c^(-ga) * wage_eff
        % donde c = wage_eff * h + r*a
        fun_corner = @(h) psi * max(h, 1e-10)^(1/frisch) - wage_eff * max(wage_eff * h + r*a, 1e-10)^(-ga);
        
        % Buscar solución entre 0 y 1
        h_init_guess = max(0.01, min(0.99, ( (wage_eff/psi) * max(wage_eff*0.5+r*a, 1e-6)^(-ga) )^frisch ));
        [h_sol, ~, flag_fzero_cnr] = fzero(fun_corner, h_init_guess, options_fzero);
        
        if flag_fzero_cnr == 1 && h_sol >= 0 && h_sol <= 1
            h_opt = h_sol;
            c_opt = wage_eff * h_opt + r*a;
        else % Si fzero falla o h_sol fuera de [0,1], la solución es h=0
            h_opt = 0;
            c_opt = r*a;
        end
    end
    
    % Asegurar consumo positivo
    c_opt = max(c_opt, 1e-9);

    % Verificar Condición KKT para el sector excluido (h_other = 0)
    % ¿MRS en h_other=0 >= Salario efectivo del sector excluido?
    % psi_other * (epsilon)^(1/frisch) >= c_opt^(-ga) * other_wage_eff
    kkt_satisfied = (other_psi * (1e-10)^(1/frisch)) >= c_opt^(-ga) * other_wage_eff;
    % Simplificando (ya que el lado izquierdo es casi 0 si frisch>0):
    kkt_satisfied = (c_opt^(-ga) * other_wage_eff <= 1e-9); % El beneficio marginal de trabajar en el otro sector es despreciable

    if kkt_satisfied
        flag = true; % Esta esquina es potencialmente óptima
        if strcmp(sector, 'formal')
            util_opt = c_opt^(1-ga)/(1-ga) - psi_f*h_opt^(1+1/frisch)/(1+1/frisch);
        else
            util_opt = c_opt^(1-ga)/(1-ga) - psi_i*h_opt^(1+1/frisch)/(1+1/frisch);
        end
    else
        flag = false; % Esta esquina no es óptima
        util_opt = -Inf;
    end

end