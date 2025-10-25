# Comparación de Modelos: Evolución del Framework de Aiyagari

## Tabla Comparativa

| Característica | Modelo Base (Moll) | Con Trabajo Endógeno | Con Dos Sectores |
|---------------|-------------------|---------------------|------------------|
| **Variables de Elección del Hogar** | | | |
| Consumo | c | c | c |
| Ahorro | s = a' - a | s = a' - a | s = a' - a |
| Trabajo Formal | - | ℓ (agregado) | ℓᶠ |
| Trabajo Informal | - | - | ℓⁱ |
| **Total Variables** | 2 | 3 | 4 |
| | | | |
| **Función de Utilidad** | | | |
| Forma | u(c) = c^(1-γ)/(1-γ) | u(c,ℓ) = c^(1-γ)/(1-γ) - ℓ^(1+1/φ)/(1+1/φ) | u(c,ℓᶠ,ℓⁱ) con desutilidades separables |
| Parámetros | γ | γ, φ | γ, φ, ψᶠ, ψⁱ |
| | | | |
| **Ingresos del Hogar** | | | |
| Salario | w·z (exógeno) | w·z·ℓ | wᶠ·z·ℓᶠ(1-τ) + wⁱ·z·ℓⁱ |
| Capital | r·a | r·a | r·a |
| Impuestos | - | - | τ sobre formal |
| | | | |
| **Firmas** | | | |
| Sectores | 1 | 1 | 2 (Formal + Informal) |
| Producción | Y = AK^α L^(1-α) | Y = AK^α L^(1-α) | Yᶠ = Aᶠ(Kᶠ)^αᶠ(Lᶠ)^(1-αᶠ), Yⁱ = Aⁱ(Lⁱ)^αⁱ |
| | | | |
| **Equilibrio** | | | |
| Mercados | K, L (exógeno) | K, L (endógeno) | K, Lᶠ, Lⁱ |
| Precios | r, w | r, w | r, wᶠ, wⁱ |
| | | | |
| **Complejidad Computacional** | | | |
| FOCs | 1 (consumo) | 2 (consumo, trabajo) | 3 (consumo, 2 trabajos) |
| Estado Estacionario | Directo | Requiere resolver ℓ₀(a,z) | Requiere punto fijo en wⁱ |
| Tiempo Ejecución | ~30 seg | ~1-2 min | ~3-5 min |

## Ecuaciones Clave de Cada Modelo

### 1. Modelo Base (Aiyagari-Moll)

**HJB:**
```
ρvⱼ(a) = max_c {u(c) + v'ⱼ(a)[wzⱼ + ra - c] + λⱼ[v₋ⱼ(a) - vⱼ(a)]}
```

**FOC:**
```
c = (v'ⱼ)^(-1/γ)
```

### 2. Modelo con Trabajo Endógeno

**HJB:**
```
ρvⱼ(a) = max_{c,ℓ} {u(c,ℓ) + v'ⱼ(a)[wzⱼℓ + ra - c] + λⱼ[v₋ⱼ(a) - vⱼ(a)]}
```

**FOCs:**
```
c = (v'ⱼ)^(-1/γ)
ℓ = (v'ⱼ·w·zⱼ)^φ
```

**Condición Implícita (ahorro cero):**
```
ℓ - (wzℓ + ra)^(-γφ)·(wz)^φ = 0
```

### 3. Modelo con Dos Sectores

**HJB:**
```
ρvⱼ(a) = max_{c,ℓᶠ,ℓⁱ} {u(c,ℓᶠ,ℓⁱ) + v'ⱼ(a)[wᶠzⱼℓᶠ(1-τ) + wⁱzⱼℓⁱ + ra - c] + λⱼ[v₋ⱼ(a) - vⱼ(a)]}
```

**FOCs:**
```
c = (v'ⱼ)^(-1/γ)
ℓᶠ = [(v'ⱼ·wᶠ·zⱼ·(1-τ))/ψᶠ]^φ
ℓⁱ = [(v'ⱼ·wⁱ·zⱼ)/ψⁱ]^φ
```

**Sistema en Estado Estacionario (s=0):**
```
c = wᶠzℓᶠ(1-τ) + wⁱzℓⁱ + ra
ψᶠ·(ℓᶠ)^(1/φ) = c^(-γ)·wᶠz(1-τ)
ψⁱ·(ℓⁱ)^(1/φ) = c^(-γ)·wⁱz
```

## Ventajas y Desventajas de Cada Modelo

### Modelo Base
**Ventajas:**
- Simple y rápido
- Bien entendido teóricamente
- Fácil calibración

**Desventajas:**
- Trabajo exógeno poco realista
- No captura respuesta laboral a políticas

### Modelo con Trabajo Endógeno
**Ventajas:**
- Captura margen intensivo laboral
- Respuestas más realistas a shocks
- Permite análisis de políticas laborales

**Desventajas:**
- Mayor complejidad computacional
- Requiere calibrar elasticidad de Frisch

### Modelo con Dos Sectores
**Ventajas:**
- Captura informalidad endógena
- Permite análisis de política fiscal
- Rica heterogeneidad en comportamiento

**Desventajas:**
- Alta complejidad computacional
- Muchos parámetros a calibrar
- Posibles equilibrios múltiples

## Aplicaciones Recomendadas

| Modelo | Mejor para estudiar |
|--------|-------------------|
| **Base** | • Distribución de riqueza<br>• Ahorro precautorio<br>• Efectos de r sobre desigualdad |
| **Trabajo Endógeno** | • Respuesta laboral a impuestos<br>• Elasticidad de oferta laboral<br>• Interacción ahorro-trabajo |
| **Dos Sectores** | • Economía informal<br>• Evasión fiscal endógena<br>• Desarrollo económico<br>• Política fiscal en economías duales |

## Parámetros Adicionales por Modelo

### Parámetros Nuevos en Trabajo Endógeno
- `φ` (frisch): Elasticidad de Frisch [0.5 - 2.0]

### Parámetros Nuevos en Dos Sectores
- `τ`: Tasa impositiva formal [0.1 - 0.4]
- `ψᶠ`: Desutilidad trabajo formal [0.8 - 1.5]
- `ψⁱ`: Desutilidad trabajo informal [0.5 - 1.5]
- `Aᶠ`: Productividad formal [calibrar]
- `Aⁱ`: Productividad informal [calibrar]
- `αᶠ`: Share capital formal [0.3 - 0.4]
- `αⁱ`: Retornos informal [0.6 - 0.9]

## Código: Transición Entre Modelos

### De Base a Trabajo Endógeno
```matlab
% Cambios principales:
% 1. Añadir elasticidad de Frisch
frisch = 0.5;

% 2. Resolver para ℓ₀ antes del loop HJB
for i=1:I
    params = [a(i), z_j, w, r, ga, frisch];
    l0(i,j) = fzero(@lab_solve_aiyagari, x0, options);
end

% 3. Modificar FOCs en upwind
lf = (dVf.*w.*zz).^frisch;
lb = (dVb.*w.*zz).^frisch;
```

### De Trabajo Endógeno a Dos Sectores
```matlab
% Cambios principales:
% 1. Separar trabajo en dos sectores
lf_f = ((dVf.*wf.*zz*(1-tau))/psi_f).^frisch;
li_f = ((dVf.*wi.*zz)/psi_i).^frisch;

% 2. Loop de punto fijo para salario informal
for iter_wi = 1:max_iter_wi
    wi = al_i * Ai * Li_old^(al_i-1);
    % ... resolver HJB y FP ...
    Li_new = integral(g * li);
end

% 3. Dos mercados de trabajo
Lf = integral(g * lf);
Li = integral(g * li);
```

## Diagnóstico de Problemas Comunes

| Síntoma | Causa Probable | Solución |
|---------|---------------|----------|
| No converge HJB | Δ muy pequeño | Aumentar Delta a 1000+ |
| Distribución explota | Grilla insuficiente | Aumentar I o amax |
| Horas > 1 | FOCs mal especificadas | Verificar restricción ℓᶠ+ℓⁱ ≤ 1 |
| Salario informal oscila | Inestabilidad punto fijo | Usar relajación: wⁱ_new = θ·wⁱ_calc + (1-θ)·wⁱ_old |
| r no converge | Exceso oferta/demanda grande | Ajustar rmin, rmax |

## Recomendaciones de Uso

1. **Para investigación sobre desigualdad básica**: Usar modelo base
2. **Para política laboral/fiscal simple**: Usar trabajo endógeno
3. **Para economías en desarrollo o informalidad**: Usar dos sectores
4. **Para calibración inicial**: Comenzar con modelo base, luego añadir complejidad
5. **Para teaching**: Progresión natural de simple a complejo

## Referencias de Implementación

- **Código Base**: https://benjaminmoll.com/codes/
- **Documentación HJB**: HACT_Numerical_Appendix.pdf
- **Trabajo Endógeno**: labor_supply.pdf (Moll)
- **Dos Sectores**: Restrepo-Echavarria (2014), Horvath (2018)
