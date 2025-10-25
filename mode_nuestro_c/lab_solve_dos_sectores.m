function [lf, li, c] = lab_solve_dos_sectores(a, z, wf, wi, r, ga, frisch, psi_f, psi_i, tau)
% LAB_SOLVE_DOS_SECTORES - Resuelve oferta laboral óptima con ahorro cero (s=0)
% para modelo con sectores formal e informal. Maneja restricciones.

% Opciones para el solver
options_fsolve = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-9, 'TolX', 1e-9);
options_fzero = optimset('Display','off', 'TolFun', 1e-9, 'TolX', 1e-9);

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
    ratio_target = (psi_i * wf * (1-tau)) / (psi_f * wi);
    
    % Función para encontrar lf en la frontera
    fun_boundary = @(lf_bound) (max(lf_bound, 1e-10) / max(1 - lf_bound, 1e-10))^(1/frisch) - ratio_target;
    
    % Resolver para lf_boundary (entre 0 y 1)
    if fun_boundary(1e-9) * fun_boundary(1-1e-9) < 0 % Si hay cambio de signo
        [lf_bound, ~, flag_fzero_bnd] = fzero(fun_boundary, [1e-9, 1-1e-9], options_fzero);
        if flag_fzero_bnd ~= 1
            lf_bound = 0.5; % Fallback si fzero falla
        end
    elseif fun_boundary(1e-9) > 0 % Solución es lf=0
        lf_bound = 0;
    else % Solución es lf=1
        lf_bound = 1;
    end
        
    lf_bound = max(0, min(1, lf_bound));
    li_bound = 1 - lf_bound;
    c_bound = wf*z*lf_bound*(1-tau) + wi*z*li_bound + r*a;
    c_bound = max(1e-9, c_bound);
    
    % Calcular utilidad en la frontera
    util_bound = c_bound^(1-ga)/(1-ga) - psi_f*lf_bound^(1+1/frisch)/(1+1/frisch) ...
                 - psi_i*li_bound^(1+1/frisch)/(1+1/frisch);
             
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