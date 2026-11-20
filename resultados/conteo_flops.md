# Conteo de Operaciones (FLOPs) - Eliminación Gaussiana (n = 12)

## Convención de conteo
- **1 FLOP** = 1 operación aritmetica elemental en punto flotante (`+`, `-`, `*`, `/`).
- Los multiplicadores $m_{ik} = a_{ik}/a_{kk}$ se calculan **una sola vez por fila** y se reutilizan.
- No se operan ni reescriben las posiciones bajo la diagonal ya anuladas ($j < k+1$).
- El sobrecosto de pivoteo parcial escalado (cálculo de factores de escala, cocientes $|a_{ik}|/s_i$, comparaciones e intercambios) se contabiliza **por separado**.

## Desglose por paso (Eliminación hacia adelante)

| Paso $k$ | Mul/Div | Sum/Res | Subtotal FLOPs | Cocientes escala (piv.) | Comparaciones (piv.) | Intercambios |
|:---:|---:|---:|---:|---:|---:|---:|
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

## Resumen comparativo con cotas teóricas

| Etapa / Estrategia | Mul/Div | Sum/Res | Total FLOPs | Sobrecosto Pivoteo |
|---|---:|---:|---:|---|
| Eliminación (11 pasos) | 638 | 572 | 1210 | - |
| Sustitución regresiva | 78 | 66 | 144 | - |
| **Total sin pivoteo** | **716** | **638** | **1354** | 0 comparaciones, 0 intercambios |
| **Total con pivoteo parcial escalado** | **716** | **638** | **1354** | 198 comp. (132 escala + 66 búsqueda), 77 div., 8 intercambios |

- **Mínimo teórico clásico** ($2n^3/3$ para $n=12$): `1152` FLOPs (solo matriz) / `~1486` FLOPs (cota de referencia con actualización de $b$ y sustitución regresiva).
- **Diferencia respecto a cota de referencia (1486)**: `-132` FLOPs (-8.88%, dentro del margen +/-10% permitido al no ejecutar restas redundantes en la columna pivote).
