# Informe técnico: eliminación gaussiana con pivoteo

## Identificación

- **Estudiantes:** Ignacio Urrea y Renato Ortiz
- **Asignatura:** Computación Numérica (INF-321)
- **Carrera:** Ingeniería Civil Informática
- **Fecha:** 20 de noviembre de 2026
- **Repositorio:** https://github.com/renaa-p/evaluacion3-pivote

## Planteamiento

Se resuelve `A x = b` con `n = 12`. La matriz tiene 144 elementos no nulos,
con magnitudes entre `1e-8` y `1e5`; por tanto, no es posible ahorrar operaciones
saltando ceros. El dato `a(1,1) = 1.3407e-08` es doce órdenes menor que
`a(2,1)`, y la fila 7 está casi en dependencia lineal de las filas 3 y 5.
La segunda característica explica el condicionamiento elevado y la primera hace
inestable comenzar sin intercambio de filas.

El objetivo numérico tiene dos componentes que deben cumplirse simultáneamente:
llevar `A` a forma triangular superior sin hacer operaciones redundantes y
conservar la mayor exactitud posible en la mantisa. El costo de la estabilidad
se controla cuantitativamente mediante el pivoteo; la cancelación catastrófica
se limita mediante el intercambio de filas y el uso de multiplicadores pequeños
y reutilizables.

Para cada fila se define la escala `s(i) = max_j |a(i,j)|`. El pivoteo
escalado compara `|a(i,k)|/s(i)`, de modo que una fila con entradas grandes no
sea elegida solo por su magnitud absoluta y las escalas incompatibles no dominen
la selección del pivote.

Los valores de `A` y `b` se conservan completos en `datos/`, sin truncar cifras.

## Estrategia de eliminación

Se ejecutan dos corridas del mismo algoritmo: sin pivoteo y con pivoteo parcial
escalado. La eliminación avanza en el orden `k = 1, 2, ..., 11`; en cada paso
se calcula una sola vez cada multiplicador `m(i,k) = a(i,k)/a(k,k)`, se guarda
en `L` y se reutiliza para actualizar las columnas `k+1..n` y el término
independiente. La parte ya anulada bajo la diagonal no se vuelve a tocar. Al
finalizar se realiza sustitución regresiva; no se calcula la inversa.

Esta elección alcanza el mínimo de trabajo de una matriz densa dentro de la
convención declarada: cada multiplicador se calcula una vez y cada posición
activa se actualiza una vez. Por eso el contador obtiene 1354 FLOPs, un 8.88 %
menos que la cota de referencia de 1486 al omitir la resta redundante de la
columna pivote; el resultado queda dentro del margen permitido de +/-10 %.

El pivote escalado maximiza `|a(i,k)|/s(i)`, donde `s(i)` es el máximo absoluto
original de la fila. El pivote inicial `a(1,1) = 1.3407e-08` se rechaza
explícitamente: aunque `a(2,1) = 1.19229739e04`, usarlo como pivote obliga a
calcular `m(2,1) = a(2,1)/a(1,1) = 8.8930960692e11`. Esa razón separa los
exponentes en `log2(|m(2,1)|) = 39.69` bits; al actualizar se amplifica el
redondeo y se pueden perder aproximadamente 40 de los 53 bits de mantisa en las
restas subsiguientes. Por eso el pivote pequeño es numéricamente inaceptable,
aunque no sea exactamente cero: el intercambio por una fila de escala compatible
evita esa división y reduce el crecimiento.

## Conteo de operaciones

Se cuenta una FLOP por cada `+`, `-`, `*` o `/` en punto flotante. El
multiplicador se cuenta una vez por fila. La asignación exacta de cero bajo la
diagonal no se cuenta como operación aritmética. Las divisiones y comparaciones
de la búsqueda de pivote se reportan aparte.

| Paso | Mul/Div | Sum/Res | Subtotal | Cocientes escala | Comparaciones | Intercambios |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 143 | 132 | 275 | 12 | 11 | 1 |
| 2 | 120 | 110 | 230 | 11 | 10 | 0 |
| 3 | 99 | 90 | 189 | 10 | 9 | 1 |
| 4 | 80 | 72 | 152 | 9 | 8 | 1 |
| 5 | 63 | 56 | 119 | 8 | 7 | 1 |
| 6 | 48 | 42 | 90 | 7 | 6 | 1 |
| 7 | 35 | 30 | 65 | 6 | 5 | 1 |
| 8 | 24 | 20 | 44 | 5 | 4 | 1 |
| 9 | 15 | 12 | 27 | 4 | 3 | 1 |
| 10 | 8 | 6 | 14 | 3 | 2 | 0 |
| 11 | 3 | 2 | 5 | 2 | 1 | 0 |
| **Eliminación** | **638** | **572** | **1210** | **77** | **66** | **8** |
| **Sustitución regresiva** | **78** | **66** | **144** | - | - | - |
| **Total** | **716** | **638** | **1354** | - | - | - |

El contador está instrumentado en `src/conteo.py` y `src/gauss.py`. El total
aritmético queda 8.88 % bajo la cota de referencia de 1486 porque no se cuenta
la resta redundante de la columna pivote, que se anula por asignación. Esto está
dentro del margen de +/-10 % indicado; el pivoteo agrega 77 divisiones, 66
comparaciones de búsqueda, 132 comparaciones de escalas iniciales y 8
intercambios, reportados por separado.

## Análisis de mantisa y pivoteo

| Estrategia | Multiplicador máximo | Factor de crecimiento | Error relativo |
|---|---:|---:|---:|
| Sin pivoteo | `2.752593e13` | `2.641224e3` | `2.380178e-10` |
| Pivoteo parcial escalado | `2.465558e9` | `1.783202` | `2.701713e-13` |

Sin pivoteo, los multiplicadores enormes fuerzan restas entre números de escalas
incompatibles: el operando pequeño pierde bits bajos al alinear exponentes. El
factor de crecimiento confirma que la eliminación amplifica los errores. El
pivoteo escalado reduce ambos efectos.

En doble precisión hay 53 bits de mantisa. La cota por condicionamiento
`kappa_2(A) = 3.084736e13` equivale a perder `log2(kappa_2(A)) = 44.81` bits
en el peor caso, dejando aproximadamente `8.19` bits efectivos, equivalentes a
`2.16` dígitos decimales. Como indicadores del crecimiento observado, los
multiplicadores máximos corresponden a `44.65` bits sin pivoteo y `31.20` bits
con pivoteo, mientras que el factor de crecimiento con pivoteo agrega solo
`log2(1.783202) = 0.83` bits. El resultado observado es mejor que la cota
porque esta es pesimista y combina peores casos que no ocurren simultáneamente.

## Condicionamiento

Con `numpy.linalg.cond` se obtiene:

- `kappa_2(A) = 3.084736e13`
- `kappa_inf(A) = 5.367827e13`
- `kappa_2(A) * epsilon = 6.849490e-03`
- Cota pesimista: aproximadamente `2.16` dígitos recuperables.

La casi dependencia lineal mencionada en el planteamiento explica el número de
condición. La cota no predice exactamente los dígitos obtenidos: es una cota
superior de peor caso. Con pivoteo, el residuo y la comparación con la solución
de alta precisión muestran entre 10 y 13 dígitos correctos antes del
refinamiento, y 16 después en esta ejecución.

## Resultados

Los dígitos afirmados provienen de comparar cada valor de doble precisión con
`mpmath` a 80 dígitos. La columna refinada se obtiene con dos pasadas de
refinamiento reutilizando `L` y `U`.

| Incógnita | Pivoteo escalado | Dígitos | Refinado | Dígitos |
|---|---:|---:|---:|---:|
| x1 | 0.7384915626999 | 13 | 0.7384915627 | 16 |
| x2 | -1.2059374168 | 12 | -1.2059374168 | 16 |
| x3 | 0.0043871952 | 10 | 0.0043871952 | 16 |
| x4 | -0.0917264538 | 11 | -0.0917264538 | 16 |
| x5 | 2.8461093725 | 12 | 2.8461093725 | 16 |
| x6 | -0.006281540899 | 10 | -0.0062815409 | 16 |
| x7 | 1.4738205916 | 12 | 1.4738205916 | 16 |
| x8 | -0.0384927164 | 11 | -0.0384927164 | 16 |
| x9 | 0.9152648037 | 12 | 0.9152648037 | 16 |
| x10 | -2.6071394852 | 12 | -2.6071394852 | 16 |
| x11 | 0.0518736294 | 11 | 0.0518736294 | 16 |
| x12 | -1.0849371625 | 12 | -1.0849371625 | 16 |

## Verificación

El residuo `r = b - A x` se calcula con `mpmath` a 80 dígitos, usando los valores
decimales originales. Las normas infinito son:

| Corrida | Norma infinito del residuo |
|---|---:|
| Sin pivoteo | `1.7717520031e-07` |
| Pivoteo escalado | `2.2908764385e-10` |
| Refinamiento (2 pasadas) | `2.6031873696e-10` |

El refinamiento mejora el error relativo de `2.701713e-13` a `4.855364e-17` y
conserva los 16 dígitos reportados; la norma del residuo no baja
monotónicamente porque la corrección se resuelve en doble precisión y el sistema
está mal condicionado. Como control independiente, `numpy.linalg.solve` produce
los mismos valores mostrados en `resultados/solucion.csv`; no se usa para la
factorización propia.

## Conclusiones

El compromiso elegido es pagar el pequeño sobrecosto del pivoteo escalado para
evitar el crecimiento de la matriz y conservar al menos 10 dígitos en las 12
incógnitas. La eliminación respeta el costo denso mínimo dentro del margen
solicitado y no hace trabajo redundante. Para `n = 1000` mantendría pivoteo,
pero usaría almacenamiento y operaciones por bloques, BLAS/LAPACK para las
actualizaciones y consideraría factorización iterativa o refinamiento en
precisión mixta.

## Reproducibilidad

- Python: `3.13.7` (compatible con Python 3.10 o superior).
- Dependencias: `numpy==2.1.3`, `mpmath==1.3.0`.
- Instalación: `pip install -r requirements.txt`.
- Ejecución: `python src/main.py`.
- Archivos reproducidos: todos los archivos de `resultados/`.
- Fuentes: enunciado de la evaluación (sección 3) y documentación oficial de
  [NumPy](https://numpy.org/doc/) y [mpmath](https://mpmath.org/doc/current/).
  Las bibliotecas se usaron para condicionamiento, control y precisión
  extendida; la factorización propia no se delega a una rutina de biblioteca.

Se utilizó un asistente de IA para apoyar la organización de la documentación,
la verificación de requisitos y la revisión de consistencia.
