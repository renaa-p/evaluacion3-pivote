from dataclasses import dataclass
from typing import List, Tuple
from src.conteo import RegistroOperaciones


@dataclass
class ResultadoGauss:
    x: List[float]
    U: List[List[float]]
    L: List[List[float]]
    perm: List[int]
    max_multiplicador: float
    factor_crecimiento: float
    max_inicial: float
    max_durante: float
    registro: RegistroOperaciones
    trazas: List[str]


def _formatear_matriz(M: List[List[float]], b_vec: List[float], titulo: str) -> str:
    lineas = [f"=== {titulo} ==="]
    for i, fila in enumerate(M):
        vals = " ".join(f"{v:14.6e}" for v in fila)
        lineas.append(f"F{i+1:02d}: {vals} | {b_vec[i]:14.6e}")
    lineas.append("")
    return "\n".join(lineas)


def sustitucion_regresiva(
    U: List[List[float]],
    b: List[float],
    registro: RegistroOperaciones | None = None,
) -> List[float]:
    n = len(U)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        suma = b[i]
        for j in range(i + 1, n):
            suma -= U[i][j] * x[j]
            if registro is not None:
                registro.sust_mul_div += 1
                registro.sust_sum_res += 1
        x[i] = suma / U[i][i]
        if registro is not None:
            registro.sust_mul_div += 1
    return x


def resolver_lu_permutado(
    L: List[List[float]],
    U: List[List[float]],
    perm: List[int],
    b_orig: List[float],
) -> List[float]:
    """Reutiliza los factores L y U sin refactorizar para resolver A * x = b_orig."""
    n = len(U)
    y = [b_orig[perm[i]] for i in range(n)]
    for k in range(n - 1):
        for i in range(k + 1, n):
            y[i] -= L[i][k] * y[k]
    return sustitucion_regresiva(U, y, registro=None)


def eliminacion_gaussiana(
    A_in: List[List[float]],
    b_in: List[float],
    pivoteo_escalado: bool = True,
) -> ResultadoGauss:
    n = len(A_in)
    A = [fila[:] for fila in A_in]
    b = b_in[:]
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    perm = list(range(n))

    estrategia = "Con pivoteo parcial escalado" if pivoteo_escalado else "Sin pivoteo"
    reg = RegistroOperaciones(estrategia=estrategia)

    max_inicial = max(abs(A[i][j]) for i in range(n) for j in range(n))
    max_durante = max_inicial
    max_mult = 0.0
    trazas: List[str] = []

    escalas = [0.0] * n
    if pivoteo_escalado:
        for i in range(n):
            mejor = abs(A[i][0])
            for j in range(1, n):
                reg.escala_inicial_comp += 1
                val = abs(A[i][j])
                if val > mejor:
                    mejor = val
            escalas[i] = mejor

    for k in range(n - 1):
        paso_reg = reg.iniciar_paso(k + 1)

        if pivoteo_escalado:
            idx_piv = k
            mejor_ratio = abs(A[k][k]) / escalas[k]
            paso_reg.piv_div += 1
            for i in range(k + 1, n):
                ratio = abs(A[i][k]) / escalas[i]
                paso_reg.piv_div += 1
                paso_reg.piv_comp += 1
                if ratio > mejor_ratio:
                    mejor_ratio = ratio
                    idx_piv = i

            if idx_piv != k:
                A[k], A[idx_piv] = A[idx_piv], A[k]
                b[k], b[idx_piv] = b[idx_piv], b[k]
                escalas[k], escalas[idx_piv] = escalas[idx_piv], escalas[k]
                perm[k], perm[idx_piv] = perm[idx_piv], perm[k]
                if k > 0:
                    L[k][:k], L[idx_piv][:k] = L[idx_piv][:k], L[k][:k]
                paso_reg.intercambios += 1

        pivote = A[k][k]
        if pivote == 0.0:
            raise ZeroDivisionError(f"Pivote nulo en el paso {k + 1}")

        for i in range(k + 1, n):
            m_ik = A[i][k] / pivote
            paso_reg.mul_div += 1
            L[i][k] = m_ik
            if abs(m_ik) > max_mult:
                max_mult = abs(m_ik)

            A[i][k] = 0.0  # Anulación exacta bajo la diagonal sin operación redundante
            for j in range(k + 1, n):
                A[i][j] -= m_ik * A[k][j]
                paso_reg.mul_div += 1
                paso_reg.sum_res += 1
                if abs(A[i][j]) > max_durante:
                    max_durante = abs(A[i][j])

            b[i] -= m_ik * b[k]
            paso_reg.mul_div += 1
            paso_reg.sum_res += 1

        trazas.append(_formatear_matriz(A, b, f"{estrategia} - Paso {k + 1}"))

    x = sustitucion_regresiva(A, b, registro=reg)
    factor_crecimiento = max_durante / max_inicial

    return ResultadoGauss(
        x=x,
        U=A,
        L=L,
        perm=perm,
        max_multiplicador=max_mult,
        factor_crecimiento=factor_crecimiento,
        max_inicial=max_inicial,
        max_durante=max_durante,
        registro=reg,
        trazas=trazas,
    )