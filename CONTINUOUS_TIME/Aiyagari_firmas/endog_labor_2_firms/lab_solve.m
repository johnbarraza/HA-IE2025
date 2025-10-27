function eq = lab_solve_dual(l,params)
% l = [lf; li] = trabajo formal e informal
% params = [a, z, wF, wI, r, gamma, frisch]

a     = params(1);
z     = params(2);
wF    = params(3);
wI    = params(4);
r     = params(5);
gamma = params(6);
phi   = params(7);  % Frisch elasticity

lf = l(1);  % trabajo formal
li = l(2);  % trabajo informal

% ingreso total (consumo)
c = wF*z*lf + wI*z*li + r*a;

% condiciones de primer orden para cada sector
eq1 = (lf)^(1/phi) * c^gamma - wF*z;
eq2 = (li)^(1/phi) * c^gamma - wI*z;

eq = [eq1; eq2];
