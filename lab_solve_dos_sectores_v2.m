function [lf, li, c] = lab_solve_dos_sectores(a, z, wf, wi, r, ga, frisch, psi_f, psi_i, tau)
% LAB_SOLVE_DOS_SECTORES - Resuelve oferta laboral óptima con ahorro cero (s=0)
% para modelo con sectores formal e informal
%
% Sistema a resolver:
% [1] c = wf*z*lf*(1-tau) + wi*z*li + r*a     (restricción presupuestaria)
% [2] psi_f * lf^(1/frisch) = c^(-ga) * wf*z*(1-tau)   (FOC trabajo formal)
% [3] psi_i * li^(1/frisch) = c^(-ga) * wi*z   (FOC trabajo informal)

% Opciones para solvers
options_fsolve = optimset('Display', 'off', 'TolFun', 1e-9, 'TolX', 1e-9, 'MaxFunEvals', 500);
options_fzero = optimset('Display', 'off', 'TolFun', 1e-9, 'TolX', 1e-9);

% Salarios efectivos
wf_eff = wf * z * (1-tau);
wi_eff = wi * z;

% =========================================================================
% PASO 1: Intentar solución interior (sin restricción lf + li <= 1)
% =========================================================================

% Función objetivo: sistema de FOCs
fun = @(x) system_focs_interior(x, a, z, wf, wi, r, ga, frisch, psi_f, psi_i, tau);

% Conjetura inicial inteligente
income_guess = r*a + 0.3*wf_eff + 0.3*wi_eff;
c_init = max(0.01, income_guess);

% De las FOCs, aproximación inicial
lf_init = min(0.5, max(0, (c_init^(-ga) * wf_eff / psi_f)^frisch));
li_init = min(0.5, max(0, (c_init^(-ga) * wi_eff / psi_i)^frisch));

x0 = [lf_init; li_init; c_init];

% Resolver
[x_sol, ~, exitflag] = fsolve(fun, x0, options_fsolve);

lf_interior = x_sol(1);
li_interior = x_sol(2);
c_interior = x_sol(3);

% Verificar si la solución es válida
valid_interior = (exitflag > 0) && ...
                 (c_interior > 0) && ...
                 (lf_interior >= -1e-9) && ...
                 (li_interior >= -1e-9) && ...
                 (lf_interior + li_interior <= 1 + 1e-9);

if valid_interior
    % Solución interior es válida
    lf = max(0, lf_interior);
    li = max(0, li_interior);
    
    % Asegurar restricción por errores numéricos
    if lf + li > 1
        total = lf + li;
        lf = lf / total;
        li = li / total;
    end
    
    c = wf*z*lf*(1-tau) + wi*z*li + r*a;
    return;
end

% =========================================================================
% PASO 2: Si solución interior falla, evaluar casos restringidos
% =========================================================================

solutions = [];

% --- Caso A: Frontera lf + li = 1 ---
if wf_eff > 0 && wi_eff > 0
    % Ratio óptimo de las FOCs: (lf/li)^(1/frisch) = (psi_i/psi_f) * (wf_eff/wi_eff)
    ratio = ((psi_i/psi_f) * (wf_eff/wi_eff))^frisch;
    
    % Con lf + li = 1: lf = ratio*li, entonces li = 1/(1+ratio)
    li_bound = 1/(1 + ratio);
    lf_bound = 1 - li_bound;
    
    % Asegurar que están en [0,1]
    lf_bound = max(0, min(1, lf_bound));
    li_bound = max(0, min(1, li_bound));
    
    c_bound = wf*z*lf_bound*(1-tau) + wi*z*li_bound + r*a;
    
    if c_bound > 0
        util_bound = utility(c_bound, lf_bound, li_bound, ga, frisch, psi_f, psi_i);
        solutions = [solutions; lf_bound, li_bound, c_bound, util_bound];
    end
end

% --- Caso B: Solo trabajo formal (li = 0) ---
if wf_eff > 0
    % FOC: psi_f * lf^(1/frisch) = c^(-ga) * wf_eff
    % donde c = wf_eff * lf + r*a
    
    fun_formal = @(lf) psi_f * lf^(1/frisch) - (wf_eff*lf + r*a)^(-ga) * wf_eff;
    
    % Buscar solución
    lf_init_corner = min(0.8, max(0.01, ((wf_eff/psi_f) * max(r*a, 0.01)^(-ga))^frisch));
    
    try
        [lf_only, ~, flag] = fzero(fun_formal, lf_init_corner, options_fzero);
        if flag > 0 && lf_only >= 0 && lf_only <= 1
            c_formal = wf_eff * lf_only + r*a;
            
            % Verificar que no es rentable trabajar en informal
            marginal_benefit_informal = c_formal^(-ga) * wi_eff;
            marginal_cost_informal = psi_i * (1e-10)^(1/frisch); % ~0
            
            if marginal_benefit_informal <= marginal_cost_informal * 10 % Margen de tolerancia
                util_formal = utility(c_formal, lf_only, 0, ga, frisch, psi_f, psi_i);
                solutions = [solutions; lf_only, 0, c_formal, util_formal];
            end
        end
    catch
        % Si falla, continuar
    end
end

% --- Caso C: Solo trabajo informal (lf = 0) ---
if wi_eff > 0
    % FOC: psi_i * li^(1/frisch) = c^(-ga) * wi_eff
    % donde c = wi_eff * li + r*a
    
    fun_informal = @(li) psi_i * li^(1/frisch) - (wi_eff*li + r*a)^(-ga) * wi_eff;
    
    % Buscar solución
    li_init_corner = min(0.8, max(0.01, ((wi_eff/psi_i) * max(r*a, 0.01)^(-ga))^frisch));
    
    try
        [li_only, ~, flag] = fzero(fun_informal, li_init_corner, options_fzero);
        if flag > 0 && li_only >= 0 && li_only <= 1
            c_informal = wi_eff * li_only + r*a;
            
            % Verificar que no es rentable trabajar en formal
            marginal_benefit_formal = c_informal^(-ga) * wf_eff;
            marginal_cost_formal = psi_f * (1e-10)^(1/frisch); % ~0
            
            if marginal_benefit_formal <= marginal_cost_formal * 10 % Margen de tolerancia
                util_informal = utility(c_informal, 0, li_only, ga, frisch, psi_f, psi_i);
                solutions = [solutions; 0, li_only, c_informal, util_informal];
            end
        end
    catch
        % Si falla, continuar
    end
end

% --- Caso D: No trabajar (lf = 0, li = 0) ---
c_nowork = r*a;
if c_nowork > 0
    util_nowork = utility(c_nowork, 0, 0, ga, frisch, psi_f, psi_i);
    solutions = [solutions; 0, 0, c_nowork, util_nowork];
end

% =========================================================================
% PASO 3: Elegir la mejor solución (mayor utilidad)
% =========================================================================

if isempty(solutions)
    % Caso de emergencia: asignar valores por defecto
    warning('No se encontró solución válida. Usando valores por defecto.');
    
    if wf_eff > wi_eff && wf_eff > 0
        lf = min(0.5, max(0.1, -r*a/wf_eff));
        li = 0;
    elseif wi_eff > 0
        lf = 0;
        li = min(0.5, max(0.1, -r*a/wi_eff));
    else
        lf = 0;
        li = 0;
    end
    
    c = max(0.001, wf*z*lf*(1-tau) + wi*z*li + r*a);
else
    % Elegir solución con mayor utilidad
    [~, best_idx] = max(solutions(:,4));
    lf = solutions(best_idx, 1);
    li = solutions(best_idx, 2);
    c = solutions(best_idx, 3);
end

% Asegurar restricciones finales
lf = max(0, min(1, lf));
li = max(0, min(1, li));

if lf + li > 1
    total = lf + li;
    lf = lf / total;
    li = li / total;
end

c = max(0.001, wf*z*lf*(1-tau) + wi*z*li + r*a);

end

%% =========================================================================
% FUNCIONES AUXILIARES
% =========================================================================

function F = system_focs_interior(x, a, z, wf, wi, r, ga, frisch, psi_f, psi_i, tau)
    % Sistema de FOCs para el caso interior
    lf = x(1);
    li = x(2);
    c = x(3);
    
    % Evitar valores problemáticos
    lf_eval = max(1e-10, lf);
    li_eval = max(1e-10, li);
    c_eval = max(1e-10, c);
    
    F = zeros(3,1);
    
    % Ecuación 1: Restricción presupuestaria con s=0
    F(1) = c - (wf*z*lf*(1-tau) + wi*z*li + r*a);
    
    % Ecuación 2: FOC trabajo formal
    F(2) = psi_f * lf_eval^(1/frisch) - c_eval^(-ga) * wf*z*(1-tau);
    
    % Ecuación 3: FOC trabajo informal
    F(3) = psi_i * li_eval^(1/frisch) - c_eval^(-ga) * wi*z;
    
    % Penalización suave para guiar la convergencia
    if lf < 0
        F(2) = F(2) + 100*lf^2;
    end
    if li < 0
        F(3) = F(3) + 100*li^2;
    end
    if c < 0
        F(1) = F(1) + 100*c^2;
    end
end

function u = utility(c, lf, li, ga, frisch, psi_f, psi_i)
    % Calcula la utilidad
    if c <= 0
        u = -Inf;
        return;
    end
    
    % Utilidad del consumo
    if ga == 1
        u_c = log(c);
    else
        u_c = c^(1-ga)/(1-ga);
    end
    
    % Desutilidad del trabajo formal
    if lf > 0
        u_lf = -psi_f * lf^(1+1/frisch)/(1+1/frisch);
    else
        u_lf = 0;
    end
    
    % Desutilidad del trabajo informal
    if li > 0
        u_li = -psi_i * li^(1+1/frisch)/(1+1/frisch);
    else
        u_li = 0;
    end
    
    u = u_c + u_lf + u_li;
    
    % Verificar NaN o Inf
    if ~isfinite(u)
        u = -Inf;
    end
end
