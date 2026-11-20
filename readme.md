# Eliminación gaussiana con pivoteo escalado

Este repositorio resuelve un sistema lineal denso de 12 incógnitas mediante
eliminación gaussiana, comparando una corrida sin pivoteo con otra con pivoteo
parcial escalado.

## Autoría y asignatura

- **Autores:** Ignacio Urrea y Renato Ortiz
- **Asignatura:** Computación Numérica (INF-321)
- **Institución:** Universidad Católica del Maule

## Requisitos

- Python 3.10 o superior.
- Instalar las dependencias fijadas desde la raíz del repositorio:

```bash
pip install -r requirements.txt
```

## Ejecución

Ejecutar desde la raíz del repositorio, es decir, desde el directorio que contiene
los archivos `readme.md`, `requirements.txt`, `informe.md` y las carpetas `datos/`,
`src/` y `resultados/`:

```bash
python src/main.py
```

El programa genera:

- `resultados/solucion.csv`: incógnitas y dígitos reportados.
- `resultados/traza_eliminacion.txt`: matrices de los 11 pasos.
- `resultados/conteo_flops.md`: conteo instrumentado de operaciones.
- `resultados/residuo.txt`: residuos, normas y métricas numéricas.

## Estructura de carpetas

```text
datos/       Matriz A y vector b de entrada.
src/         Implementación de Gauss, conteo, verificación y punto de entrada.
resultados/  Archivos generados por la ejecución.
```

## Resultado resumido

| Incógnita | Valor con pivoteo (10 dígitos) |
|---|---:|
| x1 | 0.7384915627 |
| x2 | -1.205937417 |
| x3 | 0.004387195200 |
| x4 | -0.09172645380 |
| x5 | 2.846109373 |
| x6 | -0.006281540899 |
| x7 | 1.473820592 |
| x8 | -0.03849271640 |
| x9 | 0.9152648037 |
| x10 | -2.607139485 |
| x11 | 0.05187362940 |
| x12 | -1.084937163 |

| Métrica | Resultado |
|---|---:|
| Total de FLOPs aritméticos | 1354 |
| Norma infinito del residuo con pivoteo | `2.2908764385e-10` |

## Informe técnico

El desarrollo y análisis completo están en [informe.md](informe.md).

## Nota de uso académico

Este repositorio se entrega para fines académicos. El código y los resultados
deben ser comprendidos y defendidos por el estudiante.
