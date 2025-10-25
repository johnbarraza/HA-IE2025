%Optimized for speed by SeHyoun Ahn

clear all; clc; close all;

tic; % Inicia el cronómetro

% -------------------------------------------------------------------------
% 1. DEFINICIÓN DE PARÁMETROS DEL MODELO
% -------------------------------------------------------------------------

% --- Parámetros de los Hogares ---
ga = 2;            % Coeficiente de aversión relativa al riesgo (\gamma)
rho = 0.05;        % Tasa de descuento subjetiva
d = 0.05;          % Tasa de depreciación del capital (\delta)

frisch = 1/2;     % Elasticidad de Frisch de la oferta laboral (\varphi)

% --- Parámetros de las Firmas (Producción) ---
al = 1/3;          % Participación del capital (\alpha) en Cobb-Douglas Y = A*K^al*L^(1-al)
Aprod = 0.1;       % Productividad total de los factores (A)

% --- Parámetros del Proceso Estocástico (Ingreso) ---
z1 = 1;            % Productividad baja
z2 = 2*z1;         % Productividad alta
z = [z1,z2];       % Vector de estados de productividad
la1 = 1/3;         % Tasa de Poisson z1 -> z2 (\lambda_1)
la2 = 1/3;         % Tasa de Poisson z2 -> z1 (\lambda_2)
la = [la1,la2];


tau = 0.30;       % Tasa impositiva sobre ingreso formal (ejemplo, ajustar para Perú)
nu = 1.0;         % Peso relativo del ocio/desutilidad del trabajo (ejemplo, calibrar)
theta = 0.7;      % Ratio salario informal / salario formal (wI = theta*wF) (ejemplo, calibrar)
v_inf = 0.85;     % Parámetro de rendimientos decrecientes informal (si usas producción informal)
p_audit = 0.1;    % Probabilidad de auditoría (si usas Restrepo)
phi_penalty = 1.5;% Factor de multa (si usas Restrepo)


% Trabajo agregado (L), se asume igual a la media estacionaria de 'z'
z_ave = (z1*la2 + z2*la1)/(la1 + la2);

% -------------------------------------------------------------------------
% 2. CREACIÓN DE LA GRILLA (REJILLA) DE ACTIVOS
% -------------------------------------------------------------------------

I= 1000;           % Número de puntos en la grilla de activos 'a'
amin = 0;          % Límite inferior (Restricción de no-endeudamiento, \underline{a}=0)
amax = 20;         % Límite superior de activos
% Comentario del código original: Sugerencia para reescalar 'amax' si 'Aprod' cambia
% amax = 20*(Aprod/0.1)^(1/(1-al));

a = linspace(amin,amax,I)'; % Vector columna de la grilla de activos
da = (amax-amin)/(I-1);  % Distancia entre puntos de la grilla (\Delta a)

% Matrices auxiliares para cálculos vectorizados
aa = [a,a];        % Repite la grilla de activos para ambos estados de 'z'
zz = ones(I,1)*z;  % Matriz (I x 2) con los valores de 'z'

% -------------------------------------------------------------------------
% 3. PARÁMETROS DE ITERACIÓN NUMÉRICA
% -------------------------------------------------------------------------

% --- Parámetros HJB (Bucle Interno) ---
maxit= 100;        % Número máximo de iteraciones para la HJB
crit = 10^(-6);    % Criterio de convergencia (tolerancia) para la HJB
Delta = 100;      % Tamaño del paso para el método implícito

% Inicialización de matrices
dVf = zeros(I,2);  % Derivada HJB hacia adelante
dVb = zeros(I,2);  % Derivada HJB hacia atrás
c = zeros(I,2);    % Política de consumo

% Matriz (sparse) que captura los saltos entre estados de productividad (z1 <-> z2)
Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

% --- Parámetros Equilibrio (Bucle Externo) ---
Ir = 40;           % Máximas iteraciones para encontrar 'r' (Equilibrio General)
crit_S = 10^(-5);  % Criterio de convergencia para el mercado de capitales (S=0)

% (Valores iniciales para 'r' y 'w' - serán sobreescritos)
rmax = 0.049;
r = 0.04;
w = 0.05;

% --- Parámetros de Bisección (Búsqueda de 'r') ---
r0 = 0.03;         % (No se usa, pero está definido)
rmin = 0.01;       % Límite inferior para la búsqueda de 'r'
rmax = 0.99*rho;   % Límite superior para la búsqueda de 'r' (r < rho)

% -------------------------------------------------------------------------
% 4. BUCLE EXTERNO: BÚSQUEDA DEL EQUILIBRIO GENERAL (BISECCIÓN)
% -------------------------------------------------------------------------
% Este bucle busca la tasa de interés 'r' que equilibra el mercado de capitales,
% es decir, donde la Oferta de Capital (KS) = Demanda de Capital (KD).

for ir=1:Ir;

% Almacena el historial de la búsqueda
r_r(ir)=r;
rmin_r(ir)=rmin;
rmax_r(ir)=rmax;

% --- 4.1. Resolver Precios de la Firma (Dado 'r') ---
% A partir de las FOC de la firma:
% FOC Capital: r = F_K - d = al*Aprod*(L/K)^(1-al) - d
% Despejando K (que es KD, la Demanda de Capital):
KD(ir) = (al*Aprod/(r + d))^(1/(1-al))*z_ave;
% FOC Trabajo: w = F_L = (1-al)*Aprod*(K/L)^al
w = (1-al)*Aprod*KD(ir).^al*z_ave^(-al); % Salario de equilibrio





%%% nuevo %%%%
%%% Pre-cálculo para Horas y Consumo con Ahorro Cero %%%%
options = optimset('Display','off', 'TolFun', 1e-8); % Opciones para fzero
hf0 = zeros(I,2); % Inicializa matriz para horas formales (ahorro cero)
hi0 = zeros(I,2); % Inicializa matriz para horas informales (ahorro cero)
c0 = zeros(I,2);  % Inicializa matriz para consumo (ahorro cero)
l0_total = zeros(I,2); % Horas totales (para conjetura inicial de fzero)

tic; % Inicia cronómetro
for j = 1:2 % Itera sobre los estados de productividad z1, z2
    z_val = z(j); % Productividad actual
    
    % Calcula ingresos marginales para este z
    IncomeMargF_z = w * z_val * (1 - tau);
    IncomeMargI_z = theta * w * z_val;
    
    % Determina qué sector se elige en cada punto 'a' si el ahorro es cero
    if IncomeMargF_z >= IncomeMargI_z
        sector_elegido = 'Formal';
        w_eff_z = IncomeMargF_z;
    else
        sector_elegido = 'Informal';
        w_eff_z = IncomeMargI_z;
    end
    
    % Conjetura inicial razonable para las horas (ej: basado en FOC simple)
    % c_guess ~ w_eff*0.5 + r*a (asumiendo h=0.5)
    % h_guess ~ ( (w_eff/nu) * (w_eff*0.5+r*a)^(-ga) )^frisch
    h_guess_vec = ( (w_eff_z/nu) .* max(w_eff_z*0.5 + r.*a, 1e-6).^(-ga) ).^frisch;
    h_guess_vec = max(0, min(h_guess_vec, 1)); % Limita entre 0 y 1

    for i = 1:I % Itera sobre la grilla de riqueza 'a'
        params = [a(i), z_val, w, r, ga, frisch, tau, theta, nu];
        
        % Conjetura inicial para fzero en este punto (a,z)
        h_initial_guess = max(h_guess_vec(i), 1e-4); % Evita guess de 0 exacto

        % Llama a fzero para encontrar la hora óptima en el sector elegido
        [h_opt_zero_saving, fval, exitflag] = fzero(@(h) lab_solve_h0(h, params), h_initial_guess, options);
        
        % Asegura que las horas estén entre 0 y 1
        h_opt_zero_saving = max(0, min(h_opt_zero_saving, 1));

        % Asigna las horas al sector correcto
        if strcmp(sector_elegido, 'Formal')
            hf0(i,j) = h_opt_zero_saving;
            hi0(i,j) = 0;
        else
            hf0(i,j) = 0;
            hi0(i,j) = h_opt_zero_saving;
        end
        
        % Calcula el consumo consistente con ahorro cero
        NetIncome0_ij = IncomeMargF_z * hf0(i,j) + IncomeMargI_z * hi0(i,j);
        c0(i,j) = NetIncome0_ij + r * a(i);
        
        % Guarda las horas totales (opcional, podría usarse para lmin/lmax)
        l0_total(i,j) = hf0(i,j) + hi0(i,j);

    end % fin bucle i (riqueza)
end % fin bucle j (productividad)
toc % Fin del cronómetro

lmin = l0_total(1,:); % Horas totales en amin (para condición de contorno V_a)
lmax = l0_total(I,:); % Horas totales en amax (para condición de contorno V_a)
%%% Fin Pre-cálculo %%%%


% Chequeo de que la restricción de endeudamiento es válida (natural)
if w*z(1) + r*amin < 0
    disp('CAREFUL: borrowing constraint too loose')
end

% --- 4.2. Conjetura Inicial HJB (Warm Start) ---
% Para la primera iteración (ir=1), se usa una fórmula analítica
% Para ir>1, se usa la función de valor (V) convergida de la iteración anterior
if ir==1
    v0(:,1) = (w*z(1) + r.*a).^(1-ga)/(1-ga)/rho;
    v0(:,2) = (w*z(2) + r.*a).^(1-ga)/(1-ga)/rho;
else
    v0 = V_r(:,:,ir-1); % "Warm start"
end

v = v0; % 'v' es la función de valor actual que se actualizará

% -------------------------------------------------------------------------
% 5. BUCLE INTERNO: RESOLVER HJB (Dado 'r' y 'w')
% -------------------------------------------------------------------------
for n=1:maxit
    V = v; % Guarda el valor de la iteración anterior (v^n)
    V_n(:,:,n)=V; % Almacena historial (opcional)
    
    % --- 5.1. Aproximación de derivadas (Diferencias Finitas) ---
    % Derivada hacia adelante: v'(a_i) \approx (v_{i+1} - v_i) / da
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (w*z + r.*amax).^(-ga); % Condición de contorno en amax
    
    % Derivada hacia atrás: v'(a_i) \approx (v_i - v_{i-1}) / da
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (w*z + r.*amin).^(-ga); % Condición de contorno en amin
    




% Inicializa matrices para políticas por sector y aproximación
% --- 5.2. ESQUEMA "UPWIND" con elección Formal/Informal ---

% Inicializa matrices
    cf = zeros(I,2); ssf = zeros(I,2); hf_f = zeros(I,2); hi_f = zeros(I,2);
    cb = zeros(I,2); ssb = zeros(I,2); hf_b = zeros(I,2); hi_b = zeros(I,2);
    u_f = zeros(I,2); u_b = zeros(I,2); % Utilidades para comparar

    for k = 1:2 % Itera Forward (k=1) / Backward (k=2)
        if k == 1
            dV = dVf;
        else
            dV = dVb;
        end

        % --- Calcular Óptimos Asumiendo Especialización ---

        % FOC Consumo (común a ambos casos): c = dV^(-1/ga)
        c_opt = max(dV, 1e-10).^(-1/ga);

        % --- Caso 1: Solo Formal ---
        IncomeMargF = w .* zz .* (1 - tau);
        % FOC Formal: nu * hf^(1/frisch) / c_opt^(-ga) = dV * IncomeMargF
        hf_formal_only = (dV .* IncomeMargF .* c_opt.^ga / nu).^frisch;
        hf_formal_only = max(0, min(hf_formal_only, 1));
        NetIncome_Formal = IncomeMargF .* hf_formal_only;
        % Hamiltoniano Formal: H = u(c,hf,0) + V_a*(NetIncome + ra - c)
        Util_Formal = c_opt.^(1-ga)/(1-ga) - nu.*hf_formal_only.^(1+1/frisch)/(1+1/frisch);
        Hamiltonian_F = Util_Formal + dV .* (NetIncome_Formal + r.*aa - c_opt);

        % --- Caso 2: Solo Informal ---
        % Necesitamos resolver la FOC no lineal para h_i numéricamente
        hi_informal_only = zeros(I,2); % Inicializar
        % Conjetura inicial razonable para h_i (podría mejorarse)
        h_guess_inf = ( (dV .* (1-p_audit*phi_penalty*tau).*zz.*v_inf .* c_opt.^ga / nu).^frisch ).^(1/(1+frisch*(1-v_inf)));
        h_guess_inf = max(0, min(h_guess_inf, 1));

        for j_z = 1:2 % Resuelve por columna (estado z)
            for i_a = 1:I % Resuelve por fila (estado a)
                 hi_informal_only(i_a, j_z) = solve_informal_h(c_opt(i_a, j_z), dV(i_a, j_z), z(j_z), w, r, a(i_a), ga, frisch, tau, nu, v_inf, p_audit, phi_penalty, h_guess_inf(i_a, j_z));
            end
        end
        hi_informal_only = max(0, min(hi_informal_only, 1)); % Asegura [0, 1]

        NetIncome_Informal = (1 - p_audit * phi_penalty * tau) .* zz .* hi_informal_only.^v_inf;
        % Hamiltoniano Informal: H = u(c,0,hi) + V_a*(NetIncome + ra - c)
        Util_Informal = c_opt.^(1-ga)/(1-ga) - nu.*hi_informal_only.^(1+1/frisch)/(1+1/frisch);
        Hamiltonian_I = Util_Informal + dV .* (NetIncome_Informal + r.*aa - c_opt);

        % --- Decisión del Sector (Comparar Hamiltonianos) ---
        I_choose_formal = (Hamiltonian_F >= Hamiltonian_I);
        I_choose_informal = 1 - I_choose_formal;

        % Políticas óptimas combinadas
        hf_opt = hf_formal_only .* I_choose_formal; % hf=0 si elige informal
        hi_opt = hi_informal_only .* I_choose_informal; % hi=0 si elige formal
        c_opt_final = c_opt; % El consumo es el mismo en ambos casos dado dV

        % Utilidad óptima
        Util_opt = Util_Formal .* I_choose_formal + Util_Informal .* I_choose_informal;

        % Almacenar resultados según Forward/Backward
        if k == 1
            cf = c_opt_final;
            hf_f = hf_opt;
            hi_f = hi_opt;
            NetIncome_f = NetIncome_Formal .* I_choose_formal + NetIncome_Informal .* I_choose_informal;
            ssf = NetIncome_f + r.*aa - cf; % Ahorro Forward
            u_f = Util_opt; % Utilidad Forward (para Hamiltoniano final)
        else
            cb = c_opt_final;
            hf_b = hf_opt;
            hi_b = hi_opt;
            NetIncome_b = NetIncome_Formal .* I_choose_formal + NetIncome_Informal .* I_choose_informal;
            ssb = NetIncome_b + r.*aa - cb; % Ahorro Backward
            u_b = Util_opt; % Utilidad Backward (para Hamiltoniano final)
        end
    end % Fin bucle k

    % --- Caso Ahorro Cero (I0) ---
    % ¡¡¡ ESTA PARTE TAMBIÉN NECESITA SER RECALCULADA CON PRODUCCIÓN INFORMAL !!!
    % El bloque '%%% Pre-cálculo %%%' debe adaptarse para usar la producción informal
    % y comparar utilidades (no solo Ingreso Marginal) para decidir el sector.
    % La función lab_solve_h0 debe modificarse para usar y^I = z*h^v_inf.
    % Por ahora, MANTENEMOS LA APROXIMACIÓN ANTERIOR (puede ser imprecisa):
    IncomeMargF0 = w .* zz .* (1 - tau);
    IncomeMargI0 = theta * w .* zz; % ¡Usa theta, inconsistente con producción!
    I_choose_formal0 = (IncomeMargF0 >= IncomeMargI0);
    I_choose_informal0 = 1-I_choose_formal0;
    hf0 = l0_total .* I_choose_formal0; % Usa l0_total precalculado (basado en theta)
    hi0 = l0_total .* I_choose_informal0; % Usa l0_total precalculado (basado en theta)
    NetIncome0 = IncomeMargF0 .* hf0 + IncomeMargI0 .* hi0;
    c0 = NetIncome0 + r.*aa;
    Util0 = c0.^(1-ga)/(1-ga) - nu.*(hf0+hi0).^(1+1/frisch)/(1+1/frisch);


    % --- Selección Upwind Final ---
    If = ssf > 0;
    Ib = ssb < 0;
    I0 = (1-If-Ib);

    c = cf.*If + cb.*Ib + c0.*I0; % Consumo final
    hf = hf_f.*If + hf_b.*Ib + hf0.*I0; % Horas formales finales
    hi = hi_f.*If + hi_b.*Ib + hi0.*I0; % Horas informales finales

    % --- Utilidad Final (usada en la actualización de V) ---
    u = u_f.*If + u_b.*Ib + Util0.*I0; % Usa la utilidad óptima calculada





    % --- 5.3. CONSTRUCCIÓN DE LA MATRIZ DE TRANSICIÓN 'A' ---
    % 'A' es el generador infinitesimal discretizado
    X = -min(ssb,0)/da; % Coeficiente para v_{i-1,j}
    Y = -max(ssf,0)/da + min(ssb,0)/da; % Coeficiente para v_{i,j}
    Z = max(ssf,0)/da; % Coeficiente para v_{i+1,j}
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);

    % Matriz 'A' completa (incluye saltos de 'z')
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    % Comprobación: Las filas de 'A' (el generador) deben sumar 0
    if max(abs(sum(A,2)))>10^(-9)
       disp('Improper Transition Matrix')
       %break
    end
    
    % --- 5.4. RESOLVER EL SISTEMA IMPLÍCITO ---
    % La HJB discretizada es: (1/Delta + rho)v^{n+1} - A v^{n+1} = u^n + v^n/Delta
    %                     o: B * v^{n+1} = b
    
    B = (1/Delta + rho)*speye(2*I) - A; % Matriz B
    
    % Apila (stack) los vectores para resolver el sistema lineal
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)]; % v^n
    
    b = u_stacked + V_stacked/Delta; % Vector b
    
    V_stacked = B\b; % Resuelve para v^{n+1}
    
    % Des-apila (unstack) el vector de resultado
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    % --- 5.5. Comprobación de Convergencia HJB ---
    Vchange = V - v;
    v = V; % Actualiza la función de valor

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end % Fin del bucle interno (HJB)
toc; % Fin del cronómetro para ESTA iteración de HJB

% -------------------------------------------------------------------------
% 6. FOKKER-PLANCK (KF) (Dado 'A' convergido para este 'r')
% -------------------------------------------------------------------------
% Ahora se resuelve A^T * g = 0 para encontrar la distribución estacionaria

AT = A'; % El operador de KF es la transpuesta (adjunta) del operador de HJB
b = zeros(2*I,1); % El lado derecho es cero

% "Dirty fix": El sistema A^T * g = 0 es singular (g=0 es solución).
% Se impone una condición (g_1 = 0.1) para encontrar una solución no trivial
% que luego se normaliza para que la masa total sea 1.
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row; % Reemplaza la primera fila de A^T

% Resuelve el sistema lineal modificado A^T * g = b
gg = AT\b;
% Normaliza la distribución 'gg' para que la masa total (integral) sea 1
g_sum = gg'*ones(2*I,1)*da; % Integral = suma de (g_i * da)
gg = gg./g_sum;

% Des-apila (unstack) el vector de distribución
g = [gg(1:I),gg(I+1:2*I)];

% (Opcional) Comprueba que la masa de cada tipo suma a sus probabilidades estacionarias
check1 = g(:,1)'*ones(I,1)*da;
check2 = g(:,2)'*ones(I,1)*da;

% --- 6.1. Almacenar resultados para esta iteración de 'r' ---
g_r(:,:,ir) = g; % Almacena la distribución
adot(:,:,ir) = w*zz + r.*aa - c; % Almacena la política de ahorro
V_r(:,:,ir) = V; % Almacena la función de valor

% -------------------------------------------------------------------------
% 7. ACTUALIZACIÓN DEL BUCLE DE EQUILIBRIO (BISECCIÓN)
% -------------------------------------------------------------------------

% Calcular la Oferta de Capital de los hogares K_S(r)
% K_S = \int a * g(a,z) da dz
KS(ir) = g(:,1)'*a*da + g(:,2)'*a*da; 
% Calcular el Exceso de Oferta de Capital: S(r) = K_S(r) - K_D(r)
S(ir) = KS(ir) - KD(ir); % (KD se calculó en el paso 4.1)

% --- Lógica de Bisección ---
if S(ir)>crit_S % Exceso de Oferta de capital (hogares ahorran demasiado)
    disp('Excess Supply')
    rmax = r; % El 'r' de equilibrio debe ser más bajo (para desincentivar el ahorro)
    r = 0.5*(r+rmin); % Prueba un 'r' a medio camino entre r y rmin
elseif S(ir)<-crit_S; % Exceso de Demanda de capital (hogares ahorran muy poco)
    disp('Excess Demand')
    rmin = r; % El 'r' de equilibrio debe ser más alto (para incentivar el ahorro)
    r = 0.5*(r+rmax); % Prueba un 'r' a medio camino entre r y rmax
elseif abs(S(ir))<crit_S; % ¡Equilibrio encontrado! (Exceso de Oferta/Demanda es ~0)
    display('Equilibrium Found, Interest rate =')
    disp(r)
    break % Termina el bucle externo 'ir'
end

end % Fin del bucle externo (sobre 'r')

% -------------------------------------------------------------------------
% 8. GRÁFICOS DE RESULTADOS DEL EQUILIBRIO
% -------------------------------------------------------------------------
% (El código restante es para crear y formatear los gráficos finales
% usando los resultados de la última iteración 'ir')

amax1 = 5;
amin1 = amin-0.1;

figure(1)
h1 = plot(a,adot(:,1,ir),'b',a,adot(:,2,ir),'r',linspace(amin1,amax1,I),zeros(1,I),'k--','LineWidth',2)
legend(h1,'s_1(a)','s_2(a)','Location','NorthEast')
text(-0.155,-.105,'$\underline{a}$','FontSize',16,'interpreter','latex')
line([amin amin], [-.1 .08],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Savings, $s_i(a)$','interpreter','latex')
xlim([amin1 amax1])
%ylim([-0.03 0.05])
set(gca,'FontSize',16)

figure(2)
h1 = plot(a,g_r(:,1,ir),'b',a,g_r(:,2,ir),'r','LineWidth',2)
legend(h1,'g_1(a)','g_2(a)')
text(-0.155,-.12,'$\underline{a}$','FontSize',16,'interpreter','latex')
line([amin amin], [0 max(max(g_r(:,:,ir)))],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Densities, $g_i(a)$','interpreter','latex')
xlim([amin1 amax1])
%ylim([0 0.5])
set(gca,'FontSize',16)

% Gráfico de la política de oferta laboral (l_j(a))

figure(3) % Crea una nueva figura (o usa un número diferente si ya tienes 3)
h_labor = plot(a, hf(:,1)+hi(:,1), 'b', a, hf(:,2)+hi(:,2), 'r', 'LineWidth', 2); % Grafica hf+hilegend(h_labor, 'l_1(a)', 'l_2(a)', 'Location', 'SouthWest') % Añade leyenda
grid on % Añade rejilla
xlabel('Wealth, $a$','interpreter','latex') % Etiqueta eje X
ylabel('Labor Supply, $l_i(a)$','interpreter','latex') % Etiqueta eje Y
xlim([amin1 amax1]) % Usa los mismos límites que los otros gráficos si quieres
% ylim([0.8 1.8]) % Descomenta y ajusta si necesitas cambiar los límites del eje Y
set(gca,'FontSize',16) % Ajusta tamaño de fuente

%  GRÁFICO DE LA POLÍTICA DE CONSUMO DE EQUILIBRIO
figure(4) % Crea una nueva figura (o usa un número diferente)
h_consumo = plot(a, c(:,1), 'b', a, c(:,2), 'r', 'LineWidth', 2); % Grafica c_1(a) en azul y c_2(a) en rojo
legend(h_consumo, 'c_1(a)', 'c_2(a)', 'Location', 'SouthEast') % Añade leyenda (ajusta la ubicación si es necesario)
grid on % Añade rejilla
xlabel('Wealth, $a$','interpreter','latex') % Etiqueta eje X
ylabel('Consumption, $c_i(a)$','interpreter','latex') % Etiqueta eje Y
xlim([amin1 amax1]) % Usa los mismos límites que los otros gráficos
% ylim([0 1.5]) % Descomenta y ajusta si necesitas cambiar los límites del eje Y
set(gca,'FontSize',16) % Ajusta tamaño de fuente

% Gráfico Opcional: Horas Formales
figure(5)
plot(a, hf(:,1), 'b', a, hf(:,2), 'r', 'LineWidth', 2);
legend('h_f_1(a)', 'h_f_2(a)', 'Location', 'SouthWest')
grid on; xlabel('Wealth, $a$','interpreter','latex'); ylabel('Formal Labor, $h_{f,i}(a)$','interpreter','latex'); xlim([amin1 amax1]); set(gca,'FontSize',16)

% Gráfico Opcional: Horas Informales
figure(6)
plot(a, hi(:,1), 'b', a, hi(:,2), 'r', 'LineWidth', 2);
legend('h_i_1(a)', 'h_i_2(a)', 'Location', 'NorthEast') % Puede necesitar ajustar Location
grid on; xlabel('Wealth, $a$','interpreter','latex'); ylabel('Informal Labor, $h_{i,i}(a)$','interpreter','latex'); xlim([amin1 amax1]); set(gca,'FontSize',16)