# Reporte: Analisis de Senales de Calidad de Energia mediante la Transformada de Hilbert-Huang

---

## 1. Introduccion

El analisis de senales de calidad de energia (Power Quality, PQ) es fundamental para la deteccion, clasificacion e identificacion de perturbaciones en sistemas electricos. Las senales de PQ presentan caracteristicas no estacionarias y no lineales, tales como hundimientos de tension (sags), elevaciones (swells), transitorios oscilatorios, flicker y distorsion armonica, que ocurren de manera localizada en el tiempo.

Los metodos convencionales de analisis espectral, como la Transformada Rapida de Fourier (FFT), asumen estacionariedad de la senal, lo cual limita su capacidad para capturar eventos temporales localizados. La Transformada de Hilbert-Huang (HHT), propuesta por Norden E. Huang en 1998, ofrece una alternativa adaptativa y local que no requiere funciones base predefinidas, permitiendo analizar senales no estacionarias y no lineales de manera natural.

El presente reporte tiene como objetivo evaluar el desempeno de la HHT en el analisis de senales sinteticas de calidad de energia, comparandola con la FFT y la descomposicion por paquetes Wavelet (WPT). Se analiza la capacidad de reconstruccion de cada metodo y se estudia el efecto de la ventana de analisis sobre la calidad de los resultados.

---

## 2. Fundamento Teorico

### 2.1 Transformada de Hilbert-Huang (HHT)

La HHT consta de dos etapas:

#### 2.1.1 Descomposicion Empirica en Modos (EMD)

La EMD descompone adaptativamente una senal x(t) en un conjunto finito de funciones de modo intrinseco (IMFs) y un residuo:

```
x(t) = SUM[k=1 a n] IMF_k(t) + r(t)
```

Cada IMF satisface dos condiciones:
1. El numero de extremos y el numero de cruces por cero difieren a lo sumo en uno.
2. En cualquier punto, el valor medio de la envolvente superior (maximos) y la envolvente inferior (minimos) es cero.

El proceso de extraccion de IMFs se denomina "sifting" y consiste en:
1. Identificar los maximos y minimos locales de la senal.
2. Construir las envolventes superior e inferior mediante interpolacion cubica (spline).
3. Calcular la media de las envolventes.
4. Restar la media de la senal original para obtener un candidato a IMF.
5. Repetir hasta cumplir un criterio de parada.

#### 2.1.2 Analisis Espectral de Hilbert (HSA)

Una vez obtenidas las IMFs, se aplica la Transformada de Hilbert a cada una:

```
H[IMF_k(t)] = (1/pi) * P.V. INTEGRAL[-inf a +inf] IMF_k(tau) / (t - tau) dtau
```

Esto permite obtener la senal analitica:

```
z_k(t) = IMF_k(t) + j*H[IMF_k(t)] = a_k(t) * exp(j*theta_k(t))
```

De donde se extraen:
- **Amplitud instantanea:** a_k(t) = |z_k(t)|
- **Frecuencia instantanea:** f_k(t) = (1/2*pi) * d(theta_k)/dt

#### 2.1.3 Espectro de Hilbert

El espectro de Hilbert H(omega, t) es una representacion tiempo-frecuencia que muestra la distribucion de energia:

```
H(omega, t) = SUM[k=1 a n] a_k(t)^2 * delta(omega - omega_k(t))
```

El espectro marginal h(omega) integra la energia en el tiempo:

```
h(omega) = INTEGRAL[0 a T] H(omega, t) dt
```

#### 2.1.4 Hilbert Normalizada

Para mejorar la estimacion de frecuencia instantanea en senales con modulacion de amplitud, se normaliza la IMF por su envolvente:

```
Co(t) = IMF(t) / E(t)
```

donde E(t) es la envolvente obtenida por interpolacion de los picos de |IMF(t)|. La senal normalizada Co(t) tiene amplitud unitaria, lo que permite una estimacion mas precisa de la fase y frecuencia instantanea.

### 2.2 Transformada Rapida de Fourier (FFT)

La FFT descompone una senal en componentes sinusoidales estacionarias:

```
X[k] = SUM[n=0 a N-1] x[n] * exp(-j*2*pi*k*n/N)
```

Proporciona amplitud y fase para cada frecuencia discreta f_k = k*fs/N. Su principal limitacion es que asume que las componentes frecuenciales son constantes en todo el intervalo de analisis.

### 2.3 Paquetes Wavelet (WPT)

La descomposicion por paquetes Wavelet proporciona una representacion tiempo-frecuencia mediante filtros pasa-banda en cascada. A diferencia de la Wavelet discreta (DWT), la WPT descompone tanto las aproximaciones como los detalles en cada nivel, proporcionando resolucion uniforme en frecuencia.

### 2.4 Comparacion Teorica de Metodos

| Caracteristica | FFT | HHT (EMD) | WPT |
|----------------|-----|-----------|-----|
| Base | Funciones predefinidas (senos/cosenos) | Adaptativa (derivada de los datos) | Predefinida (wavelet madre) |
| Estacionariedad | Requiere | No requiere | No requiere |
| Resolucion temporal | Ninguna (promedio global) | Instantanea | Variable segun nivel |
| Resolucion frecuencial | fs/N (fija) | Variable (instantanea) | Dependiente del nivel |
| Linealidad | Requiere | No requiere | Requiere |
| Reconstruccion | Perfecta (IFFT) | Perfecta (suma de IMFs) | Perfecta |

---

## 3. Metodologia

### 3.1 Senales de Prueba

Se utilizaron tres senales sinteticas generadas con el toolbox PQSyDa (Power Quality Synthetic Data), con los siguientes parametros base:

- **Frecuencia fundamental:** 60 Hz
- **Frecuencia de muestreo:** 18,000 Hz
- **Duracion:** 10 ciclos (0.1667 s)
- **Muestras totales:** 3,000

#### Senal #5: Transitorio Oscilatorio

```
x(t) = A*sin(2*pi*f*t - phi) + A*beta*exp(-(t-t1)/tau)*sin(2*pi*fn*(t-t1) - theta)*u(t)
```

Parametros generados:
- Frecuencia del transitorio: fn = 634.2 Hz
- Constante de tiempo: tau = 20.1 ms
- Amplitud del transitorio: beta = 0.22
- Intervalo del transitorio: [0.040, 0.076] s

#### Senal #18: Armonicos + Swell + Flicker

```
x(t) = A*(1+beta*u(t))*(alpha1*sin(2*pi*f*t) + alpha3*sin(3*2*pi*f*t) + alpha5*sin(5*2*pi*f*t))*(1+lambda*sin(2*pi*ff*t))
```

Parametros generados:
- Armonico 3ro: alpha3 = 13.9%
- Armonico 5to: alpha5 = 8.7%
- Flicker: lambda = 0.066, ff = 17.8 Hz
- Swell: beta = 0.41, intervalo [0.015, 0.159] s

#### Senal #28: Swell + Armonicos + Flicker + Transitorio Oscilatorio

```
x(t) = A*(sin(2*pi*f*t) + (harm*beta*u + osc_transient)*flicker)
```

Parametros generados:
- Armonicos: alpha3 = 6.2%, alpha5 = 13.4%
- Flicker: lambda = 0.089, ff = 14.9 Hz
- Swell: beta = 0.46, intervalo [0.046, 0.069] s
- Oscilatorio: fn = 606.6 Hz, tau = 20.3 ms

### 3.2 Ventanas de Analisis

Cada senal se analizo con dos ventanas:

| Ventana | Muestras | Duracion | Ciclos de 60 Hz | Rango (readmatrix) |
|---------|----------|----------|-----------------|-------------------|
| Completa | 3,000 | 0.167 s | 10 | [2 1 3001 2] |
| Reducida | 1,400 | 0.078 s | ~4.7 | [2 1 1401 2] |

### 3.3 Procesamiento Aplicado

El script de MATLAB implementa las siguientes etapas:

1. **EMD:** Descomposicion de la senal en IMFs mediante la funcion `emd()` de MATLAB.
2. **HSA:** Calculo de amplitud y frecuencia instantanea mediante la Transformada de Hilbert.
3. **Espectro de Hilbert:** Construccion de la representacion H(omega,t) y el espectro marginal h(omega).
4. **Hilbert Normalizada:** Normalizacion de la IMF1 por su envolvente para mejorar la estimacion de frecuencia.
5. **FFT:** Espectro de Fourier para comparacion.
6. **WPT:** Descomposicion por paquetes Wavelet (nivel 8, wavelet db40).
7. **Reconstruccion:** Comparacion de la reconstruccion FFT (picos detectados) vs HHT (suma de IMFs + residuo).

### 3.4 Metricas de Evaluacion

- **RMSE (Root Mean Square Error):** Cuantifica el error de reconstruccion.
- **Numero de componentes:** Cantidad de elementos necesarios para la reconstruccion.
- **Estabilidad de frecuencia instantanea:** Evaluacion cualitativa de la consistencia temporal.

---

## 4. Interpretacion de Resultados

### 4.1 Descomposicion EMD

#### Efecto de la ventana de analisis sobre las IMFs

| Senal | 3000 muestras | 1400 muestras |
|-------|---------------|---------------|
| #5 Oscillatory | 4 IMFs. Separacion clara: IMF1 = transitorio, IMF2 = fundamental 60 Hz | 4 IMFs. Separacion temporal: IMF1 = parte estacionaria, IMF2 = zona del transitorio |
| #18 HarmSwellFlicker | 5 IMFs. IMF1 = fundamental+armonicos, IMF2 = swell (amplitud creciente), Residuo = tendencia DC | 4 IMFs. Mode mixing presente — las componentes se contaminan entre si |
| #28 SwellHarmFlickerOsc | 4 IMFs. IMF1 = transitorio+fundamental, IMF2 = swell, Residuo = tendencia descendente | 4 IMFs. Fuerte mode mixing — dificultad para separar componentes |

**Hallazgo clave:** Con 3000 muestras, el EMD logra separar correctamente los eventos temporales (transitorio) de las componentes estacionarias (fundamental). Con 1400 muestras, se presenta mode mixing significativo en las senales complejas (#18 y #28).

#### Mode Mixing

El mode mixing es el principal problema observado. Ocurre cuando:
- Una misma IMF contiene oscilaciones de escalas temporales muy diferentes.
- Oscilaciones similares se distribuyen en diferentes IMFs.

Este fenomeno se acentua con:
- Ventanas de analisis cortas (menos ciclos para la interpolacion de envolventes).
- Senales con modulacion de amplitud (swell, flicker).
- Multiples perturbaciones simultaneas.

### 4.2 Espectro de Hilbert y Espectro Marginal

#### Identificacion de la fundamental (60 Hz)

| Senal | 3000 muestras | 1400 muestras |
|-------|---------------|---------------|
| #5 | Pico claro en 60 Hz (energia ~0.35) | Pico en 60 Hz (energia ~0.14) |
| #18 | Energia dominante en ~5 Hz (NO identifica 60 Hz) | Energia en ~5-15 Hz (NO identifica 60 Hz) |
| #28 | Picos distribuidos en 35-85 Hz | Energia en ~10 Hz (NO identifica 60 Hz) |

**Hallazgo clave:** El espectro marginal de Hilbert tiene dificultades con senales que presentan modulacion de amplitud (swell + flicker). La modulacion de amplitud es interpretada por el EMD como una oscilacion de baja frecuencia, desplazando la energia hacia frecuencias menores a la fundamental real.

#### Representacion tiempo-frecuencia H(omega, t)

El espectro de Hilbert muestra ventajas sobre la FFT en la localizacion temporal de eventos:
- En la senal #5, el transitorio oscilatorio aparece como una concentracion de energia en alta frecuencia (>300 Hz) limitada al intervalo [0.04, 0.07] s.
- Fuera del transitorio, se observa una linea estable en ~60 Hz.

### 4.3 Hilbert Normalizada

La normalizacion de la IMF1 por su envolvente mostro resultados variables:

| Senal | 3000 muestras | 1400 muestras |
|-------|---------------|---------------|
| #5 | Frecuencia estable ~60 Hz con spike puntual en la transicion | Cambio de frecuencia visible: 60 Hz -> 350 Hz (transitorio) |
| #18 | 60 Hz + spikes en zona del transitorio | Inestable (0-400 Hz) — no interpretable |
| #28 | ~60 Hz estable en toda la ventana | Completamente inestable |

**Hallazgo clave:** La Hilbert Normalizada funciona bien cuando la envolvente E(t) es suave y no cruza por valores cercanos a cero. Cuando E(t) tiene transiciones abruptas (como en la frontera de un transitorio), la division IMF/E(t) amplifica el ruido, generando spikes en la frecuencia instantanea.

La senal #28 con 3000 muestras mostro el resultado mas satisfactorio: la normalizacion logro extraer ~60 Hz estable a pesar de ser la senal mas compleja, porque con suficientes datos la IMF1 captura bien la fundamental con envolvente suave.

### 4.4 Reconstruccion de Senales

#### Reconstruccion HHT

| Senal | Ventana | IMFs | RMSE |
|-------|---------|------|------|
| #5 | 1400 | 4 | 6.37 x 10^-17 V |
| #5 | 3000 | 4 | 6.21 x 10^-17 V |
| #18 | 1400 | 4 | 1.46 x 10^-16 V |
| #18 | 3000 | 5 | 9.97 x 10^-17 V |
| #28 | 1400 | 4 | 1.46 x 10^-16 V |
| #28 | 3000 | 4 | 2.21 x 10^-16 V |

La reconstruccion HHT siempre produce RMSE del orden de 10^-16 a 10^-17 V (precision de punto flotante de doble precision). Esto se debe a que la EMD es una descomposicion completa: x(t) = SUM(IMFs) + Residuo es una identidad matematica, no una aproximacion.

#### Reconstruccion FFT

| Senal | Ventana | Componentes | RMSE |
|-------|---------|-------------|------|
| #5 | 3000 | 2 | 0.0343 V |
| #5 | 3000 | 2 | 0.1291 V |
| #18 | 3000 | 5 | 0.0891 V |
| #18 | 1400 | 7 | 0.3861 V |
| #28 | 1400 | 13 | 0.4652 V |
| #28 | 3000 | 2 | 0.1291 V |

**Observaciones:**
1. **Senal #5 (simple):** Solo 2 componentes logran RMSE = 0.034 V. El error se concentra exclusivamente en el intervalo del transitorio oscilatorio. Fuera de ese intervalo, la reconstruccion es casi perfecta.
2. **Senal #18 (intermedia):** Con 3000 muestras y 5 componentes (60, 180, 200, 300, 360 Hz), se logra RMSE = 0.089 V. La FFT captura correctamente los armonicos estables.
3. **Senal #28 (compleja):** Con 1400 muestras, aun usando 13 componentes, el RMSE = 0.465 V. La FFT no puede representar los eventos no estacionarios (swell, flicker, transitorio) sin importar cuantas componentes se utilicen.

**Hallazgo clave:** El RMSE de la reconstruccion FFT no depende solo del numero de componentes, sino de la naturaleza de la senal. Para senales no estacionarias, aumentar componentes tiene rendimientos decrecientes porque la FFT fundamentalmente no puede representar variaciones temporales de amplitud.

### 4.5 Comparacion FFT vs HHT vs WPT

#### Espectro FFT
- Proporciona picos puntuales y claros en las frecuencias exactas de las componentes armonicas.
- Ideal para identificar frecuencias cuando la senal es estacionaria.
- Falla con eventos transitorios (la energia se distribuye en multiples bins).

#### Espectro Marginal HHT
- Proporciona espectros mas "suaves" y anchos.
- Refleja la naturaleza no estacionaria — el ancho del pico indica variabilidad frecuencial.
- Problema: la modulacion de amplitud desplaza la energia a frecuencias artificialmente bajas.

#### Espectro WPT
- Proporciona buena resolucion en bandas de frecuencia.
- Las amplitudes son proporcionadas en bandas (no frecuencias puntuales).
- Mejor compromiso tiempo-frecuencia que la FFT para senales no estacionarias.
- Limitacion: la resolucion depende del nivel de descomposicion y la wavelet elegida.

### 4.6 Efecto de la Ventana de Analisis

| Aspecto | 3000 muestras (0.167 s) | 1400 muestras (0.078 s) |
|---------|-------------------------|-------------------------|
| EMD - Separacion de modos | Buena | Pobre (mode mixing) |
| Espectro marginal | Picos mas definidos | Picos difusos |
| Hilbert Normalizada | Frecuencia estable | Inestable |
| FFT - RMSE | 0.03 - 0.13 V | 0.34 - 0.47 V |
| FFT - Resolucion frecuencial | 6 Hz (fs/N) | 12.9 Hz (fs/N) |
| Ciclos de fundamental capturados | 10 | ~4.7 |

**Limites minimos observados para analisis confiable:**
- Senal simple (#5): 1400 muestras (~4.7 ciclos) son suficientes para EMD basico.
- Senales complejas (#18, #28): Se requieren al menos 3000 muestras (~10 ciclos) para separacion adecuada de modos.
- La Hilbert Normalizada requiere envolventes suaves, lo cual se logra mejor con mas datos.

### 4.7 Capacidad de Extrapolacion

Se evaluo la posibilidad de reconstruir la senal completa a partir de una ventana parcial:

| Metodo | Capacidad de extrapolacion | Condicion |
|--------|---------------------------|-----------|
| FFT | Si — puede extrapolar indefinidamente | Solo para componentes estacionarias (armonicos constantes) |
| HHT | No — solo reconstruye dentro de la ventana | Las IMFs son funciones locales sin modelo parametrico |

Para la senal #5 con ventana reducida: la FFT extrapola correctamente la fundamental de 60 Hz fuera de la ventana, pero no predice el transitorio oscilatorio si este queda fuera del intervalo analizado.

---

## 5. Conclusiones

### 5.1 Sobre la Transformada de Hilbert-Huang

1. **Reconstruccion perfecta:** La HHT (suma de IMFs + residuo) siempre reconstruye la senal original con error del orden de la precision numerica de la maquina (~10^-16 V). Esto es una propiedad intrinseca del metodo, no una ventaja analitica, ya que la EMD es una descomposicion completa (identidad matematica).

2. **Separacion tiempo-frecuencia:** La principal ventaja de la HHT es su capacidad para separar eventos temporales localizados (transitorios) de componentes estacionarias (fundamental, armonicos). Esto se demostro claramente en la senal #5, donde el EMD aislo el transitorio oscilatorio en una IMF separada de la fundamental.

3. **Limitaciones con modulacion de amplitud:** El EMD presenta dificultades con senales que combinan swell + flicker, ya que la modulacion de amplitud es interpretada como una oscilacion de baja frecuencia, desplazando el espectro marginal hacia frecuencias inferiores a la fundamental real.

4. **Sensibilidad a la ventana de analisis:** El EMD requiere suficientes ciclos de la componente mas lenta para funcionar correctamente. Con ventanas cortas, se presenta mode mixing que compromete la interpretacion fisica de las IMFs.

### 5.2 Sobre la Comparacion de Metodos

5. **FFT vs HHT para reconstruccion:** La comparacion de RMSE entre FFT y HHT para reconstruccion no es una comparacion justa, ya que la HHT usa todas las componentes (descomposicion completa) mientras la FFT solo usa los picos detectados. Una comparacion mas significativa seria usar solo las IMFs mas energeticas vs los picos FFT.

6. **Complementariedad de metodos:** La FFT es superior para identificar frecuencias exactas en senales estacionarias y para extrapolar senales fuera de la ventana de analisis. La HHT es superior para localizar eventos en el tiempo y analizar senales no estacionarias. La WPT ofrece un compromiso intermedio.

7. **Recomendacion practica:** Para senales de PQ con perturbaciones mixtas, se recomienda un enfoque combinado:
   - FFT para identificar las componentes armonicas estacionarias.
   - HHT para localizar temporalmente los eventos transitorios y no estacionarios.
   - WPT como verificacion cruzada de la distribucion de energia en bandas.

### 5.3 Sobre la Ventana Minima de Analisis

8. **Minimo practico observado:** Para senales con fundamental de 60 Hz y perturbaciones simples (un solo tipo de evento), 1400 muestras (~4.7 ciclos, 0.078 s) son suficientes para un analisis basico. Para senales con multiples perturbaciones simultaneas, se requieren al menos 3000 muestras (~10 ciclos, 0.167 s).

9. **Compromiso datos vs calidad:** Reducir la ventana de analisis impacta principalmente:
   - La resolucion frecuencial de la FFT (Df = fs/N).
   - La capacidad del EMD para separar modos (mode mixing).
   - La suavidad de la envolvente en la Hilbert Normalizada.

10. **Estrategias de mitigacion:** Para ventanas cortas se recomienda: extension por espejo de la senal, uso de EEMD/CEEMDAN en lugar de EMD clasico, y seleccion de wavelet mas corta para WPT.

---

## Referencias

- Huang, N.E. et al. (1998). "The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis." Proceedings of the Royal Society A, 454, 903-995.
- PQSyDa Toolbox: github.com/Micke1995/PQSyDa
- IEEE Std 1159-2019. "IEEE Recommended Practice for Monitoring Electric Power Quality."

---

*Reporte generado con base en el analisis experimental de senales sinteticas de calidad de energia mediante MATLAB R2024.*
