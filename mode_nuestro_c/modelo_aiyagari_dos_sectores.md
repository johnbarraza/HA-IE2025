# Modelo Aiyagari con Sectores Formal e Informal
## Extensión con Agentes Heterogéneos y Trabajo Endógeno

### 1. Estructura del Modelo

#### 1.1 Problema del Hogar

Los hogares heterogéneos resuelven:

```
max E₀ ∫₀^∞ e^(-ρt) u(cₜ, ℓₜᶠ, ℓₜⁱ) dt

s.t.
ȧₜ = wᶠₜzₜℓₜᶠ(1-τ) + wⁱₜzₜℓₜⁱ + rₜaₜ - cₜ
aₜ ≥ ā
ℓₜᶠ + ℓₜⁱ ≤ 1
```

Donde:
- `ℓᶠ`: horas trabajadas en sector formal
- `ℓⁱ`: horas trabajadas en sector informal
- `wᶠ`: salario formal
- `wⁱ`: salario informal
- `τ`: tasa impositiva sobre trabajo formal
- `z ∈ {z₁, z₂}`: productividad idiosincrática (proceso Poisson)

#### 1.2 Función de Utilidad

Proponemos dos opciones:

**Opción A: Consumo agregado con desutilidad separable del trabajo**
```
u(c, ℓᶠ, ℓⁱ) = c^(1-γ)/(1-γ) - ψᶠ(ℓᶠ)^(1+1/φ)/(1+1/φ) - ψⁱ(ℓⁱ)^(1+1/φ)/(1+1/φ)
```

**Opción B: Consumos diferenciados (à la Restrepo)**
```
u(cᶠ, cⁱ, ℓᶠ, ℓⁱ) = [(cᶠ)^ε + β(cⁱ)^ε]^((1-γ)/ε)/(1-γ) - ψᶠ(ℓᶠ)^(1+1/φ)/(1+1/φ) - ψⁱ(ℓⁱ)^(1+1/φ)/(1+1/φ)
```

#### 1.3 Ecuación HJB

La ecuación Hamilton-Jacobi-Bellman será:

```
ρvⱼ(a) = max_{c,ℓᶠ,ℓⁱ} {u(c,ℓᶠ,ℓⁱ) + v'ⱼ(a)[wᶠzⱼℓᶠ(1-τ) + wⁱzⱼℓⁱ + ra - c] + λⱼ[v₋ⱼ(a) - vⱼ(a)]}
```

#### 1.4 Condiciones de Primer Orden

Para la **Opción A** (consumo agregado):

1. **Consumo**: `uᶜ = v'ⱼ(a)` → `c = (v'ⱼ)^(-1/γ)`

2. **Trabajo formal**: `-uₗᶠ = v'ⱼ(a)wᶠzⱼ(1-τ)` → `ψᶠ(ℓᶠ)^(1/φ) = v'ⱼ(a)wᶠzⱼ(1-τ)`

3. **Trabajo informal**: `-uₗⁱ = v'ⱼ(a)wⁱzⱼ` → `ψⁱ(ℓⁱ)^(1/φ) = v'ⱼ(a)wⁱzⱼ`

Esto implica:
```
ℓᶠ = [(v'ⱼ(a)wᶠzⱼ(1-τ))/ψᶠ]^φ
ℓⁱ = [(v'ⱼ(a)wⁱzⱼ)/ψⁱ]^φ
```

### 2. Firmas

#### 2.1 Sector Formal
Función de producción Cobb-Douglas:
```
Yᶠ = Aᶠ(Kᶠ)^αᶠ(Lᶠ)^(1-αᶠ)
```

FOCs:
- Capital: `r + δ = αᶠAᶠ(Lᶠ/Kᶠ)^(1-αᶠ)`
- Trabajo: `wᶠ = (1-αᶠ)Aᶠ(Kᶠ/Lᶠ)^αᶠ`

#### 2.2 Sector Informal
Función de producción con retornos decrecientes (solo trabajo):
```
Yⁱ = Aⁱ(Lⁱ)^αⁱ,  αⁱ < 1
```

FOC trabajo: `wⁱ = αⁱAⁱ(Lⁱ)^(αⁱ-1)`

### 3. Equilibrio

#### 3.1 Mercados de Trabajo
```
Lᶠ = ∫∫ ℓᶠ(a,z) g(a,z) da dz
Lⁱ = ∫∫ ℓⁱ(a,z) g(a,z) da dz
```

#### 3.2 Mercado de Capital
```
K = ∫∫ a g(a,z) da dz = Kᶠ
```

#### 3.3 Mercado de Bienes
- **Opción A**: `Y = Yᶠ + pYⁱ = C + δK`
- **Opción B**: `Yᶠ = Cᶠ + δK` y `Yⁱ = Cⁱ`

### 4. Algoritmo Numérico

#### 4.1 Pasos del Algoritmo

1. **Inicialización**:
   - Conjetura inicial para `r`
   - Grilla de activos `a`
   - Estados de productividad `z = [z₁, z₂]`

2. **Loop Externo (Equilibrio General)**:
   
   a. Dado `r`, calcular precios de equilibrio:
      - `Kᶠ` desde FOC de capital formal
      - `wᶠ` desde FOC de trabajo formal
      - Resolver para `Lⁱ` y `wⁱ` (requiere punto fijo)

   b. **Resolver HJB (Loop Interno)**:
      - Para cada `(a,z)`, resolver óptimos `ℓᶠ`, `ℓⁱ`, `c`
      - Construir matriz de transición `A`
      - Método implícito para actualizar `v`

   c. **Resolver Fokker-Planck**:
      - Encontrar distribución estacionaria `g`

   d. **Verificar equilibrio de mercados**:
      - Calcular `KS = ∫ a·g da`
      - Si `|KS - Kᶠ| > tol`, actualizar `r` (bisección)

#### 4.2 Desafíos Computacionales

1. **Problema de punto fijo en salarios informales**:
   - `wⁱ` depende de `Lⁱ`
   - `Lⁱ` depende de decisiones individuales `ℓⁱ`
   - Las decisiones dependen de `wⁱ`

2. **Condiciones de ahorro cero**:
   - Resolver implícitamente para `(ℓᶠ₀, ℓⁱ₀)` tal que:
     ```
     wᶠz(1-τ)ℓᶠ₀ + wⁱzℓⁱ₀ + ra = c₀
     ```
   - Con FOCs de trabajo

### 5. Extensiones y Consideraciones

#### 5.1 Riesgo de Detección (à la Restrepo)
Añadir probabilidad `p` de ser detectado trabajando informalmente:
```
Ingreso informal esperado = (1-p)wⁱzℓⁱ - p·Multa
```

#### 5.2 Capital de Trabajo (à la Horvath)
Sector formal necesita financiar fracción `η` del salario:
```
wᶠ_efectivo = wᶠ(1 + ηr)
```

#### 5.3 Restricciones de Horas
Podríamos imponer:
- Límite total: `ℓᶠ + ℓⁱ ≤ ℓ_max < 1`
- Límites sectoriales: `ℓᶠ ≤ ℓᶠ_max`, `ℓⁱ ≤ ℓⁱ_max`

### 6. Calibración Sugerida

| Parámetro | Descripción | Valor Base |
|-----------|-------------|------------|
| γ | Aversión al riesgo | 2 |
| ρ | Tasa de descuento | 0.05 |
| φ | Elasticidad Frisch | 0.5 |
| τ | Tasa impositiva formal | 0.2-0.3 |
| αᶠ | Share capital formal | 0.33 |
| αⁱ | Retornos informal | 0.7-0.8 |
| ψᶠ/ψⁱ | Ratio desutilidad | 0.8-1.2 |
| Aᶠ, Aⁱ | Productividades | Calibrar para tamaños sectoriales objetivo |

### 7. Outputs Esperados

El modelo debería generar:
1. **Distribución de riqueza** más desigual que modelo base
2. **Asignación sectorial** endógena según productividad
3. **Informalidad** decreciente en riqueza
4. Los más pobres trabajan más en informal
5. Los más ricos trabajan principalmente en formal
6. **Movilidad** entre sectores según shocks de productividad

### 8. Validación

Comparar con datos de:
- Tamaño del sector informal (% PIB, % empleo)
- Distribución de ingresos por sector
- Transiciones laborales entre sectores
- Correlación riqueza-formalidad
