# Guía de Reproducción: Calibración ABC-SMC con $\beta$ Condicionado

Este documento detalla las instrucciones paso a paso para reproducir los experimentos de calibración con $\beta$ condicionado dinámicamente según la tasa de proliferación celular ($d$), tanto para la ventana de retención clínica ($15\text{--}35\%$) como para la banda de conservación basal ($75\text{--}100\%$) y la comparativa de los 8 casos analíticos (A a H).

---

## 1. Obtención del Código y Ramas

El código está disponible y sincronizado en el repositorio remoto:

```bash
git fetch origin
git checkout experiment/retention-75-100-victor
```

*(Alternativamente, la rama base con el ajuste de 15–35% está en `experiment/homeostasis-beta-conditioned`).*

---

## 2. Compilación de los Binarios C++

Asegúrate de compilar los ejecutables de simulación en modo `Release`:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target run_bootstrapping run_stage_trajectories run_paper_features -j$(nproc)
```

---

## 3. Estructura de Scripts Principales

Dentro del directorio `experiments/exp08_cs_abc_smc_homeostasis/`:

1. **`abc_smc.py`**:
   - Motor de calibración ABC-SMC.
   - En las líneas 50–57 implementa la derivación analítica de la banda de retención:
     $$\beta \in \left[1 - \frac{hi^{1/80}}{(1 - \mu)(1 + d)},\; 1 - \frac{lo^{1/80}}{(1 - \mu)(1 + d)}\right]$$
   - Muestrea dinámicamente $\beta$ en función de la tasa de proliferación $d$ de cada partícula y el valor de $\mu$ (por defecto $\mu = 0.003$).

2. **`run_pipeline.sh`**:
   - Script automatizado que orquesta en un solo paso:
     1. Calibración ABC-SMC (8 generaciones, 50 partículas).
     2. Simulación de trayectorias de estadios celulares con los parámetros MAP.
     3. Renderizado de las figuras de publicación (`Figure3`, `Figure_stage_dynamics`).
   - Cada ejecución se almacena en una subcarpeta versionada con timestamp dentro de `experiments/exp08_cs_abc_smc_homeostasis/runs/`.

3. **`victor_8cases/simulate_victor_8cases.py`**:
   - Simula en paralelo las 8 combinaciones (Casos A a H) propuestas por Víctor:
     - 300 réplicas Monte Carlo para evaluar la penetrancia acumulada frente a Kuchenbaecker (edades 30 a 80).
     - 100 tejidos virtuales para trazar la retención de células basales y la acumulación tumoral.
     - Genera una figura comparativa de 4 paneles en `victor_8cases/figures/Figure_victor_8cases_comparison.png`.

---

## 4. Comandos de Ejecución

### Opción A: Banda de Retención 75–100% (Propuesta de Víctor)
Impone que a los 80 años se conserve entre el 75% y el 100% de la población basal sana:

```bash
bash experiments/exp08_cs_abc_smc_homeostasis/run_pipeline.sh \
  --n-cells 600 \
  --particles 50 \
  --n-sim 100 \
  --generations 8 \
  --healthy-min 0.75 \
  --healthy-max 1.00 \
  --tag "retention_75_100"
```

> **Resultado observado:** Al retener una gran masa de células basales vivas con $\mu = 0.003$, la penetrancia a los 80 años alcanza el 92.0% (por encima del 72% clínico), con un error mínimo estancado en Weighted SSE $\approx 69.45$.

---

### Opción B: Banda de Retención 15–35% (Prueba de Degradación Progresiva)
Reproduce la dinámica observada en la reunión del jueves (degradación biológicamente suave de las células basales hasta $\sim 20\%$ a los 80 años):

```bash
bash experiments/exp08_cs_abc_smc_homeostasis/run_pipeline.sh \
  --n-cells 600 \
  --particles 50 \
  --n-sim 100 \
  --generations 8 \
  --healthy-min 0.15 \
  --healthy-max 0.35 \
  --tag "retention_15_35"
```

> **Resultado observado:** Ajuste casi perfecto a la penetrancia de Kuchenbaecker (todas las edades dentro del IC 95%, Weighted SSE $= 1.05$), con $d \approx 0.135$, $\beta \approx 0.134$.

---

### Opción C: Comparativa Exhaustiva de los 8 Casos Analíticos (A a H)
Evalúa hacia adelante los 8 escenarios sin recalibrar los parámetros apoptóticos:

```bash
python3 experiments/exp08_cs_abc_smc_homeostasis/victor_8cases/simulate_victor_8cases.py
```

Resultados y figuras se guardan en:
- `experiments/exp08_cs_abc_smc_homeostasis/victor_8cases/data/victor_8cases_penetrance_summary.csv`
- `experiments/exp08_cs_abc_smc_homeostasis/victor_8cases/figures/Figure_victor_8cases_comparison.png`
