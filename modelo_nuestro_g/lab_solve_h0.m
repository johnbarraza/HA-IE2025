function eq = lab_solve_h0(h_sector, params)
% Resuelve la hora óptima (h_sector) en el sector elegido
% cuando el ahorro es cero (c = NetIncome + r*a).

% Desempaqueta parámetros
a = params(1);
z_val = params(2); % Productividad z para este agente
w = params(3);
r = params(4);
ga = params(5);
frisch = params(6);
tau = params(7); % Tasa impositiva formal
theta = params(8);% Ratio salario informal/formal
nu = params(9);   % Peso del ocio/desutilidad trabajo

% Calcula ingresos marginales netos
IncomeMargF = w * z_val * (1 - tau);
IncomeMargI = theta * w * z_val;

% Determina el ingreso marginal efectivo (del sector elegido)
w_eff = max(IncomeMargF, IncomeMargI);

% Calcula el consumo implícito por la restricción de ahorro cero
% c = w_eff * h_sector + r*a
c_implicit = w_eff * h_sector + r * a;

% Asegura que el consumo sea positivo (evita errores si r*a es muy negativo)
if c_implicit <= 1e-10
    % Si el consumo implícito no es positivo, la FOC no se puede satisfacer
    % Retorna un valor grande para indicar que h_sector no es una buena solución.
    % Nota: Podrías necesitar un manejo más robusto si amin es muy negativo.
     eq = 1e10; 
     return;
end

% Evalúa la FOC del trabajo: MRS = V_a * IngresoMarginalEfectivo
% nu * (h_sector)^(1/frisch) / c_implicit^(-ga) = c_implicit^(-ga) * w_eff
% Simplificando: nu * (h_sector)^(1/frisch) = c_implicit^(-ga) * w_eff

% Reordena para fzero (forma eq = 0):
% h_sector^(1/frisch) - (w_eff / nu) * c_implicit^(-ga) = 0
eq = h_sector^(1/frisch) - (w_eff / nu) * c_implicit^(-ga);

% Alternativamente, la forma que dejé en el pensamiento:
% eq = h_sector - ((w_eff / nu) * c_implicit^(-ga))^frisch;
% Ambas son equivalentes. La primera puede ser numéricamente más estable.

end