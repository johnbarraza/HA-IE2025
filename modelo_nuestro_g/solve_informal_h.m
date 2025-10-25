function h_i_opt = solve_informal_h(c_val, Va_val, z_val, w_val, r_val, a_val, ga, frisch, tau, nu, v_inf, p_audit, phi_penalty, h_guess)
% Resuelve la FOC no lineal para las horas informales óptimas, h_i,
% dados el consumo (c_val) y la derivada del valor (Va_val).

% Define la función cuya raíz es la h_i óptima
% FOC: nu * hi^(1/frisch) / c_val^(-ga) = Va_val * (1 - p*phi*tau) * z * v_inf * hi^(v_inf - 1)
% Reordenando a f(h_i) = 0:
% nu * hi^(1/frisch + 1 - v_inf) - Va_val * c_val^ga * (1-p*phi*tau) * z_val * v_inf = 0

params_foc = [c_val, Va_val, z_val, ga, frisch, tau, nu, v_inf, p_audit, phi_penalty];
foc_func = @(h) nu * max(h, 1e-10).^(1/frisch + 1 - v_inf) - Va_val * c_val^ga * (1-p_audit*phi_penalty*tau) * z_val * v_inf;

% Opciones para el solver
options_foc = optimset('Display','off', 'TolFun', 1e-7);

% Resuelve para h_i
[h_i_opt, fval, exitflag] = fzero(foc_func, max(h_guess, 1e-6), options_foc);

% Asegura que esté entre 0 y 1
h_i_opt = max(0, min(h_i_opt, 1));

% Manejo de casos donde fzero podría fallar
if exitflag ~= 1 || isnan(h_i_opt)
    % Si fzero falla, podría ser que la solución óptima sea h_i = 0
    % Comprueba si MRS < Ingreso Marginal en h_i cercano a 0
    h_test = 1e-6;
    mrs_at_0 = nu * h_test^(1/frisch) / c_val^(-ga);
    income_marg_at_0 = Va_val * (1 - p_audit * phi_penalty * tau) * z_val * v_inf * h_test^(v_inf - 1);
    if mrs_at_0 < income_marg_at_0
         h_i_opt = 0; % La solución es la esquina h_i=0
    else
         warning('Fzero falló al encontrar h_i óptimo.');
         h_i_opt = 0; % O asigna un valor por defecto o detiene
    end
end

end