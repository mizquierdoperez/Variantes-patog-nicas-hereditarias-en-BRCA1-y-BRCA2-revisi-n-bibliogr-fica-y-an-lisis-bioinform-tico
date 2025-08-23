# -*- coding: utf-8 -*-
"""
Created on Sat Aug 16 10:55:42 2025

@author: Mauricio Izquierdo
"""
import pandas as pd
import requests
import time
import os
import re
import matplotlib.pyplot as plt

# ---------------- CONFIGURACIÓN ----------------
VEP_URL = "https://rest.ensembl.org/vep/human/hgvs"

# Archivo de entrada y salida - Para repetibilidad, sustituir por ruta propia.
INPUT_EXCEL = "C:/Users/Mauricio Izquierdo/Desktop/Master_Bioinformática/mutaciones_BRCA-TFM.xlsx"
OUTPUT_EXCEL = "C:/Users/Mauricio Izquierdo/Desktop/Master_Bioinformática/resultados_prediccion.xlsx"

# Información para contruir el HGVS completo
TRANSCRITOS = {
    "BRCA1": "NM_007294.3",
    "BRCA2": "NM_000059.4"
}

SLEEP_TIME = 0.5


# ---------------- FUNCIONES ----------------
def clean_hgvs_version(hgvs_string):
    """Elimina la versión del transcrito (ej. .3) del HGVS."""
    if hgvs_string:
        return re.sub(r'\.\d+', '', hgvs_string)
    return None

def query_vep(hgvs_completo):
    """
    Consulta la API de VEP utilizando el HGVS completo en el cuerpo JSON de una petición POST.
    """
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    data = {"hgvs_notations": [hgvs_completo]}
    try:
        r = requests.post(VEP_URL, headers=headers, json=data, timeout=30)
        if r.status_code == 200:
            return r.json()
        else:
            print(f"[ERROR] VEP falló {r.status_code} para {hgvs_completo}. Respuesta: {r.text}")
            return {}
    except requests.exceptions.RequestException as e:
        print(f"[ERROR] Excepción de VEP para {hgvs_completo}: {e}")
        return {}

def parse_vep_result(doc):
    """
    Parsea el resultado de la API de VEP de forma segura, buscando el transcrito canónico.
    """
    result = {}
    try:
        if isinstance(doc, list) and len(doc) > 0:
            # Encuentra el transcrito canónico si existe
            canonical_tc = next((t for t in doc[0].get("transcript_consequences", []) if t.get("canonical") == 1), None)
            
            if canonical_tc:
                tc_to_use = canonical_tc
            else:
                # Si no hay transcrito canónico, usa el primer transcrito disponible
                tc_list = doc[0].get("transcript_consequences", [])
                tc_to_use = tc_list[0] if tc_list else None

            if tc_to_use:
                result["vep_consequence"] = ";".join(tc_to_use.get("consequence_terms", []))
                result["vep_impact"] = tc_to_use.get("impact")
    except Exception as e:
        print(f"[ERROR] Fallo al parsear resultado de VEP: {e}")
    return result


def build_hgvs(gen, variante_c):
    """
    Crea el HGVS completo a partir del gen y la notación de la variante, eliminando la versión.
    """
    transcript = TRANSCRITOS.get(str(gen).strip().upper())
    if transcript:
        # Se elimina la versión del transcrito para ambas APIs
        transcript_sin_version = clean_hgvs_version(transcript)
        return f"{transcript_sin_version}:{str(variante_c).strip()}"
    else:
        print(f"[ADVERTENCIA] Gen {gen} no reconocido, variante {variante_c} no procesada")
        return None


# ---------------- PROCESO ----------------
try:
    df = pd.read_excel(INPUT_EXCEL, header=2)
except FileNotFoundError:
    print(f"[ERROR] El archivo de entrada no se encuentra en la ruta: {INPUT_EXCEL}")
    exit()

if "Gen" not in df.columns or "Variante_HGVS" not in df.columns:
    raise ValueError("El Excel debe contener columnas 'Gen' y 'Variante_HGVS'.")

resultados = []

for index, row in df.iterrows():
    gen = str(row["Gen"]).strip()
    variante_c = str(row["Variante_HGVS"]).strip()

    hgvs_completo = build_hgvs(gen, variante_c)
    
    fila = {"Gen": gen, "Variante_HGVS": variante_c, "HGVS_completo": hgvs_completo}

    if hgvs_completo:
        print(f"[{index + 1}/{len(df)}] Procesando variante: {hgvs_completo}")
        
        vep_data = query_vep(hgvs_completo)
        vep_parsed = parse_vep_result(vep_data)
        fila.update(vep_parsed)
        
    resultados.append(fila)
    print(f"  --> Resultado para {hgvs_completo}: {fila}")
    time.sleep(SLEEP_TIME)

# Guardar resultados
df_result = pd.DataFrame(resultados)
df_result.to_excel(OUTPUT_EXCEL, index=False)

print(f"Proceso completado. Resultados guardados en '{OUTPUT_EXCEL}'.")


#---------- GRÁFICAS COMPARATIVAS ------------

# Definimos una función que nos ayude a crear todas las gráficas de la misma forma mas adelante
def plot_percentage_bar(df, category_col, group_col, title):
    """
    df: DataFrame con los datos
    category_col: columna con la variable a graficar (p.ej. 'Diagnóstico asociado')
    group_col: columna para agrupar (p.ej. 'Gen')
    title: título de la gráfica
    """
    # Reemplazar NaN por "NA"
    df[category_col] = df[category_col].fillna("NA")

    # Calcular porcentajes
    counts = df.groupby([group_col, category_col]).size().reset_index(name="count")
    totals = df.groupby(group_col).size().reset_index(name="total")
    merged = counts.merge(totals, on=group_col)
    merged["percentage"] = merged["count"] / merged["total"] * 100

    # Pivot para graficar
    pivot_df = merged.pivot(index=group_col, columns=category_col, values="percentage").fillna(0)

    # Plot
    ax = pivot_df.plot(kind="bar", stacked=True, figsize=(8,6))
    plt.ylabel("Porcentaje (%)")
    plt.title(title)
    plt.legend(title=category_col, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    # Añadir etiquetas dentro de las barras
    for container in ax.containers:
        labels = [f"{w:.1f}%" if w > 0 else "" for w in container.datavalues]  # las etiquetas de menos del 0% no salen
        ax.bar_label(
            container,
           labels = labels,
            label_type="center",
            fontsize=8,
            color="white",
            weight="bold"
        )

    plt.show()

    

# Graficamos los resultados de nuestra recision bibliográfica

# Patogenicidad (Diagnóstico asociado)
plot_percentage_bar(df, "Diagnóstico asociado", "Gen",
                    "Distribución de la patogenicidad (revisión bibliográfica)")

# Tipo de mutación
plot_percentage_bar(df, "Tipo de mutación ", "Gen",
                    "Distribución del tipo de mutación (revisión bibliográfica)")

# Graficamos los resultados de nuestras consultas realizadas al VEP mediante el codigo anterior

# Patogenicidad según impacto de VEP
plot_percentage_bar(df_result, "vep_impact", "Gen",
                    "Distribución de la patogenicidad según VEP")

# Tipo de mutación según VEP consequence
plot_percentage_bar(df_result, "vep_consequence", "Gen",
                    "Distribución del tipo de mutación según VEP")


# ------------- MAPEO Y COMPARACIÓN ------------ 

# Mapeos
impact_map = {
    "HIGH": "Patogénica",
    "MODERATE": "Ligeramente patogénica",
    "MODIFIER": "Ligeramente patogénica"
}

mutation_map = {
    "frameshift_variant": "Frameshift",
    "missense_variant": "Missense",
    "splice_region_variant": "Splice",
    "splice_acceptor_variant": "Splice",
    "splice_donor_variant": "Splice",
    "inframe_insertion": "Nonsense",
    "intron_variant": "Nonsense",
    "start_lost": "Nonsense",
    "stop_gained": "Nonsense"
}


# Función de comparación
def comparar_y_graficar(df, df_result, col_var="Variante_HGVS"):
    # Unimos ambos conjuntos por la variante
    merged = df.merge(df_result, on=col_var, suffixes=("_df", "_res"))

    resultados = []
    for _, row in merged.iterrows():
        gen = row["Gen_df"]

        # Comparar Patogenicidad 
        patog_df = row["Diagnóstico asociado"]
        impact_res = str(row.get("vep_impact", ""))

        coincide_pato = False
        if pd.notna(patog_df) and impact_res:
            mapped = impact_map.get(impact_res, None)
            if mapped and mapped.lower() in patog_df.lower():
                coincide_pato = True

        # Comparar Tipo de Mutación
        tipo_df = row["Tipo de mutación "]
        conseq_res = str(row.get("vep_consequence", ""))

        coincide_tipo = False
        if pd.notna(tipo_df) and conseq_res:
            conseq_list = conseq_res.split(";")
            mapped_cons = {mutation_map[c] for c in conseq_list if c in mutation_map}
            if tipo_df in mapped_cons:
                coincide_tipo = True

        resultados.append({
            "Gen": gen,
            "Coincide_patogenicidad": coincide_pato,
            "Coincide_tipo": coincide_tipo
        })

    comp_df = pd.DataFrame(resultados)


    # Graficamos directamente con la llamada de la función

    for col, titulo in [
        ("Coincide_patogenicidad", "Concordancia en patogenicidad"),
        ("Coincide_tipo", "Concordancia en tipo de mutación")
    ]:
        counts = comp_df.groupby(["Gen", col]).size().reset_index(name="count")
        totals = comp_df.groupby("Gen").size().reset_index(name="total")
        merged_counts = counts.merge(totals, on="Gen")
        merged_counts["percentage"] = merged_counts["count"] / merged_counts["total"] * 100

        pivot_df = merged_counts.pivot(index="Gen", columns=col, values="percentage").fillna(0)

        ax = pivot_df.plot(kind="bar", stacked=True, figsize=(7,5))
        plt.ylabel("Porcentaje (%)")
        plt.title(titulo)
        plt.legend(title="Coincide", labels=["No", "Sí"], bbox_to_anchor=(1.05, 1), loc="upper left")

        # Etiquetas dentro de las barras
        for container in ax.containers:
            labels = [f"{w:.1f}%" if w > 0 else "" for w in container.datavalues]
            ax.bar_label(container, labels=labels, label_type="center", fontsize=8, color="white", weight="bold")

        plt.tight_layout()
        plt.show()

comparar_y_graficar(df, df_result, col_var="Variante_HGVS")














