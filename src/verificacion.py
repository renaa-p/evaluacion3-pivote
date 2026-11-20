from dataclasses import dataclass
from typing import List, Tuple
import math
import mpmath as mp
import numpy as np
from src.gauss import ResultadoGauss, resolver_lu_permutado


mp.mp.dps = 80  # 80 dígitos decimales de precisión extendida


@dataclass
class ReporteVerificacion:
    kappa_2: float
    kappa_inf: float
    digitos_cota_pesimista: float
    residuo_sin_piv: List[mp.mpf]
    norma_inf_sin_piv: mp.mpf
    residuo_con_piv: List[mp.mpf]
    norma_inf_con_piv: mp.mpf
    x_refinado: List[float]
    residuo_refinado: List[mp.mpf]
    norma_inf_refinado: mp.mpf
    x_numpy: List[float]
    x_exacto_mp: List[mp.mpf]
    digitos_correctos_sin_piv: List[int]
    digitos_correctos_con_piv: List[int]
    digitos_correctos_refinado: List[int]
    error_relativo_sin_piv: float
    error_relativo_con_piv: float
    error_relativo_refinado: float


def _residuo_extendido(
    A_str: List[List[str]],
    b_str: List[str],
    x_vec: List[float],
) -> Tuple[List[mp.mpf], mp.mpf]:
    n = len(b_str)
    A_mp = [[mp.mpf(val) for val in fila] for fila in A_str]
    b_mp = [mp.mpf(val) for val in b_str]
    x_mp = [mp.mpf(val) for val in x_vec]

    r = []
    for i in range(n):
        ax_i = mp.fsum(A_mp[i][j] * x_mp[j] for j in range(n))
        r.append(b_mp[i] - ax_i)
    norma_inf = max(abs(val) for val in r)
    return r, norma_inf


def _contar_digitos_significativos(aprox: float, exacto: mp.mpf) -> int:
    if exacto == 0:
        return 16 if aprox == 0 else 0
    err_rel = abs(mp.mpf(aprox) - exacto) / abs(exacto)
    if err_rel == 0:
        return 16
    digs = -float(mp.log10(err_rel))
    return max(0, min(16, int(math.floor(digs))))


def ejecutar_verificacion(
    A_str: List[List[str]],
    b_str: List[str],
    res_sin: ResultadoGauss,
    res_con: ResultadoGauss,
) -> ReporteVerificacion:
    n = len(b_str)
    A_np = np.array([[float(v) for v in fila] for fila in A_str], dtype=np.float64)
    b_np = np.array([float(v) for v in b_str], dtype=np.float64)

    # Condicionamiento en norma 2 e infinito
    kappa_2 = float(np.linalg.cond(A_np, 2))
    kappa_inf = float(np.linalg.cond(A_np, np.inf))
    eps_ieee = float(np.finfo(np.float64).eps)
    cota_error = kappa_2 * eps_ieee
    digitos_cota = max(0.0, -math.log10(cota_error))

    # Control con biblioteca estándar
    x_numpy = np.linalg.solve(A_np, b_np).tolist()

    # Solucion de referencia en 80 dígitos decimales (mpmath)
    A_mp_mat = mp.matrix([[mp.mpf(v) for v in fila] for fila in A_str])
    b_mp_vec = mp.matrix([mp.mpf(v) for v in b_str])
    x_exacto_mp = list(mp.lu_solve(A_mp_mat, b_mp_vec))

    # Residuos en precisiopn extendida
    r_sin, norm_sin = _residuo_extendido(A_str, b_str, res_sin.x)
    r_con, norm_con = _residuo_extendido(A_str, b_str, res_con.x)

    # Refinamiento iterativo (2 pasadas reutilizando L y U de res_con)
    x_ref = res_con.x[:]
    for _ in range(2):
        r_actual, _ = _residuo_extendido(A_str, b_str, x_ref)
        r_float = [float(val) for val in r_actual]
        delta_x = resolver_lu_permutado(res_con.L, res_con.U, res_con.perm, r_float)
        x_ref = [x_ref[i] + delta_x[i] for i in range(n)]

    r_ref, norm_ref = _residuo_extendido(A_str, b_str, x_ref)

    dig_sin = [_contar_digitos_significativos(res_sin.x[i], x_exacto_mp[i]) for i in range(n)]
    dig_con = [_contar_digitos_significativos(res_con.x[i], x_exacto_mp[i]) for i in range(n)]
    dig_ref = [_contar_digitos_significativos(x_ref[i], x_exacto_mp[i]) for i in range(n)]

    norm_x_exact = float(mp.sqrt(mp.fsum(v**2 for v in x_exacto_mp)))
    err_sin = float(mp.sqrt(mp.fsum((mp.mpf(res_sin.x[i]) - x_exacto_mp[i])**2 for i in range(n)))) / norm_x_exact
    err_con = float(mp.sqrt(mp.fsum((mp.mpf(res_con.x[i]) - x_exacto_mp[i])**2 for i in range(n)))) / norm_x_exact
    err_ref = float(mp.sqrt(mp.fsum((mp.mpf(x_ref[i]) - x_exacto_mp[i])**2 for i in range(n)))) / norm_x_exact

    return ReporteVerificacion(
        kappa_2=kappa_2,
        kappa_inf=kappa_inf,
        digitos_cota_pesimista=digitos_cota,
        residuo_sin_piv=r_sin,
        norma_inf_sin_piv=norm_sin,
        residuo_con_piv=r_con,
        norma_inf_con_piv=norm_con,
        x_refinado=x_ref,
        residuo_refinado=r_ref,
        norma_inf_refinado=norm_ref,
        x_numpy=x_numpy,
        x_exacto_mp=x_exacto_mp,
        digitos_correctos_sin_piv=dig_sin,
        digitos_correctos_con_piv=dig_con,
        digitos_correctos_refinado=dig_ref,
        error_relativo_sin_piv=err_sin,
        error_relativo_con_piv=err_con,
        error_relativo_refinado=err_ref,
    )