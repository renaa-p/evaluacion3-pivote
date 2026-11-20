import csv
import os
import sys

# Asegura importaciones limpias al ejecutar `python src/main.py` desde la raíz
RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if RAIZ not in sys.path:
    sys.path.insert(0, RAIZ)

from src.gauss import eliminacion_gaussiana
from src.conteo import generar_reporte_markdown
from src.verificacion import ejecutar_verificacion


def cargar_datos(ruta_a: str, ruta_b: str):
    A_str, A_flt = [], []
    with open(ruta_a, "r", encoding="utf-8") as f:
        for fila in csv.reader(f):
            if fila:
                limpia = [x.strip() for x in fila]
                A_str.append(limpia)
                A_flt.append([float(x) for x in limpia])

    b_str, b_flt = [], []
    with open(ruta_b, "r", encoding="utf-8") as f:
        for fila in csv.reader(f):
            if fila:
                val = fila[0].strip()
                b_str.append(val)
                b_flt.append(float(val))

    return A_str, A_flt, b_str, b_flt


def main() -> None:
    ruta_a = os.path.join(RAIZ, "datos", "matriz_a.csv")
    ruta_b = os.path.join(RAIZ, "datos", "vector_b.csv")
    dir_res = os.path.join(RAIZ, "resultados")
    os.makedirs(dir_res, exist_ok=True)

    A_str, A_flt, b_str, b_flt = cargar_datos(ruta_a, ruta_b)

    # 1. Corridas sin pivoteo y con pivoteo parcial escalado
    res_sin = eliminacion_gaussiana(A_flt, b_flt, pivoteo_escalado=False)
    res_con = eliminacion_gaussiana(A_flt, b_flt, pivoteo_escalado=True)

    # 2. Verificación independiente en precisión extendida y refinamiento
    ver = ejecutar_verificacion(A_str, b_str, res_sin, res_con)

    # 3. Guardar resultados/traza_eliminacion.txt
    ruta_traza = os.path.join(dir_res, "traza_eliminacion.txt")
    with open(ruta_traza, "w", encoding="utf-8") as f:
        f.write("TRAZA DE ELIMINACIÓN GAUSSIANA (11 PASOS)\n")
        f.write("=" * 60 + "\n\n")
        f.write("\n".join(res_con.trazas))
        f.write("\n" + "=" * 60 + "\n")
        f.write("\n".join(res_sin.trazas))

    # 4. Guardar resultados/conteo_flops.md
    ruta_conteo = os.path.join(dir_res, "conteo_flops.md")
    generar_reporte_markdown(res_sin.registro, res_con.registro, ruta_conteo)

    # 5. Guardar resultados/solucion.csv (dígitos defendibles en pivoteo y tras refinamiento)
    ruta_sol = os.path.join(dir_res, "solucion.csv")
    with open(ruta_sol, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "incógnita",
            "valor_pivoteo_escalado",
            "dígitos_defendibles_pivoteo",
            "valor_refinado",
            "dígitos_defendibles_refinado",
            "valor_sin_pivoteo",
            "dígitos_correctos_sin_pivoteo",
            "control_numpy",
        ])
        for i in range(12):
            d_piv = ver.digitos_correctos_con_piv[i]
            d_ref = ver.digitos_correctos_refinado[i]
            d_sin = ver.digitos_correctos_sin_piv[i]
            val_piv_fmt = f"{res_con.x[i]:.{d_piv}g}"
            val_ref_fmt = f"{ver.x_refinado[i]:.{d_ref}g}"
            val_sin_fmt = f"{res_sin.x[i]:.{d_sin}g}"
            writer.writerow([
                f"x{i+1}",
                val_piv_fmt,
                d_piv,
                val_ref_fmt,
                d_ref,
                val_sin_fmt,
                d_sin,
                f"{ver.x_numpy[i]:.12g}",
            ])

    # 6. Guardar resultados/residuo.txt
    ruta_residuo = os.path.join(dir_res, "residuo.txt")
    with open(ruta_residuo, "w", encoding="utf-8") as f:
        f.write("VERIFICACIÓN NUMÉRICA Y RESIDUO EN PRECISIÓN EXTENDIDA (80 DÍGITOS)\n")
        f.write("=" * 72 + "\n")
        f.write(f"Condicionamiento kappa_2(A)   : {ver.kappa_2:.6e}\n")
        f.write(f"Condicionamiento kappa_inf(A) : {ver.kappa_inf:.6e}\n")
        f.write(f"Cota pesimista kappa_2 * eps  : {ver.kappa_2 * 2.220446049250313e-16:.6e} "
                f"(~{ver.digitos_cota_pesimista:.2f} dígitos)\n\n")

        f.write("--- MÉTRICAS DE ESTABILIDAD ---\n")
        f.write(f"[Sin pivoteo]    Máximo multiplicador: {res_sin.max_multiplicador:.6e} | "
                f"Factor crecimiento: {res_sin.factor_crecimiento:.6e} | "
                f"Error rel. ||dx||/||x||: {ver.error_relativo_sin_piv:.6e}\n")
        f.write(f"[Con pivoteo]    Máximo multiplicador: {res_con.max_multiplicador:.6e} | "
                f"Factor crecimiento: {res_con.factor_crecimiento:.6e} | "
                f"Error rel. ||dx||/||x||: {ver.error_relativo_con_piv:.6e}\n")
        f.write(f"[Refinado (2x)]  Error rel. ||dx||/||x||: {ver.error_relativo_refinado:.6e}\n\n")

        f.write("--- NORMA INFINITO DEL RESIDUO r = b - Ax ---\n")
        f.write(f"||r||_inf (Sin pivoteo)           : {float(ver.norma_inf_sin_piv):.10e}\n")
        f.write(f"||r||_inf (Pivoteo escalado)      : {float(ver.norma_inf_con_piv):.10e}\n")
        f.write(f"||r||_inf (Tras refinamiento)     : {float(ver.norma_inf_refinado):.10e}\n\n")

        f.write("--- VECTOR RESIDUO r = b - Ax (PIVOTEO ESCALADO) ---\n")
        for i, val in enumerate(ver.residuo_con_piv):
            f.write(f"r[{i+1:02d}] = {float(val):+.12e}\n")

        f.write("\n--- VECTOR RESIDUO r = b - Ax (TRAS REFINAMIENTO ITERATIVO) ---\n")
        for i, val in enumerate(ver.residuo_refinado):
            f.write(f"r_ref[{i+1:02d}] = {float(val):+.12e}\n")

    print("Ejecución completada. Revisa la carpeta resultados/ para los valores exactos.")


if __name__ == "__main__":
    main()