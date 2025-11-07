import pandas as pd
import networkx as nx
import numpy as np
import argparse
import os
import mygene


def convert_hugo_to_entrez(gene_list):
    """
    Convierte una lista de símbolos de genes HUGO a IDs de Entrez usando MyGene.info.
    """
    mg = mygene.MyGeneInfo()
    print("Convirtiendo símbolos HUGO a IDs de Entrez...")
    # 'scopes': tipo de identificador origen
    # 'fields': tipo de identificador salida
    # 'species': para humanos
    results = mg.querymany(gene_list, scopes='symbol', fields='entrezgene', species='human')

    entrez_ids = []
    mapping = {res['query']: res.get('entrezgene') for res in results if res.get('entrezgene')}

    for gene in gene_list:
        entrez_id = mapping.get(gene)
        if entrez_id:
            # Los IDs en la red probablemente son strings, así que los convertimos a string
            entrez_ids.append(str(entrez_id))
            print(f"  '{gene}' -> '{entrez_id}'")
        else:
            print(f"  Advertencia: No se encontró un ID de Entrez para '{gene}'. Será ignorado.")

    return entrez_ids


def run_rwr(graph, seeds, alpha=0.85, max_iter=100, tolerance=1e-6):
    """
    Implementa el algoritmo de Random Walk with Restart (RWR), el motor
    principal detrás de herramientas de propagación como GUILD.

    Args:
        graph (nx.Graph): El grafo de NetworkX. Representación computacional de la red biológica, construida a partir del archivo de texto de entrada (e.g., network_guild.txt) mediante la librería NetworkX
        seeds (list): Lista de nodos semilla.
        alpha (float): Probabilidad de reinicio (teletransportación a las semillas).
        max_iter (int): Número máximo de iteraciones.
        tolerance (float): Umbral para determinar la convergencia.

    Returns:
        pd.DataFrame: Un DataFrame con los nodos y sus puntuaciones finales, ordenado.
    """

    nodes = list(graph.nodes())
    # Importante: Asegurarse de que los nodos del grafo también sean strings para la comparación
    nodes_str = [str(n) for n in nodes]
    node_to_idx = {node: i for i, node in enumerate(nodes_str)}

    adj_matrix = nx.to_numpy_array(graph, nodelist=nodes, weight='weight', dtype=float)

    with np.errstate(divide='ignore', invalid='ignore'):
        degree = adj_matrix.sum(axis=0)
        transition_matrix = adj_matrix / degree
        transition_matrix[np.isnan(transition_matrix)] = 0

    num_nodes = len(nodes)
    p0 = np.zeros(num_nodes)
    for seed in seeds:
        if seed in node_to_idx:
            p0[node_to_idx[seed]] = 1.0

    if p0.sum() == 0:
        print("Error crítico: Ninguna de las semillas convertidas existe en el grafo. Abortando RWR.")
        return None
    p0 /= p0.sum()

    pt = np.copy(p0)

    print("\nIniciando Random Walk with Restart...")
    for i in range(max_iter):
        pt_next = (1 - alpha) * np.dot(transition_matrix, pt) + alpha * p0
        diff = np.linalg.norm(pt_next - pt, 1)
        print(f"Iteración {i + 1}, Cambio: {diff:.6f}")
        if diff < tolerance:
            print(f"Convergencia alcanzada en la iteración {i + 1}.")
            pt = pt_next
            break
        pt = pt_next

    results_df = pd.DataFrame({'gene_id': nodes_str, 'score': pt})
    return results_df.sort_values(by='score', ascending=False).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description="Realiza propagación en red con conversión automática de IDs.")
    parser.add_argument("--network", required=True, help="Ruta al archivo de la red (ej. network_guild.txt).")
    parser.add_argument("--seeds", required=True, help="Ruta al archivo de genes semilla (en formato HUGO).")
    parser.add_argument("--output", required=True, help="Ruta para guardar el archivo de resultados.")
    args = parser.parse_args()

    # --- Carga de datos ---
    print(f"Cargando la red desde '{args.network}'...")
    try:
        network_df = pd.read_csv(args.network, sep=r'\s+', header=None, names=['protein1', 'protein2', 'weight'])
    except FileNotFoundError:
        print(f"Error: No se pudo encontrar el archivo de red en '{args.network}'")
        return

    print(f"Cargando genes semilla desde '{args.seeds}'...")
    try:
        seeds_df = pd.read_csv(args.seeds, header=None, names=['gene'])
        # Limpiamos espacios en blanco de las semillas HUGO
        hugo_seeds = seeds_df['gene'].str.strip().tolist()
    except FileNotFoundError:
        print(f"Error: No se pudo encontrar el archivo de semillas en '{args.seeds}'")
        return

    # --- CONVERSIÓN ---
    entrez_seeds = convert_hugo_to_entrez(hugo_seeds)
    if not entrez_seeds:
        print("Error: No se pudo convertir ninguna semilla a ID de Entrez. Abortando.")
        return

    # --- Construcción del Grafo ---
    print("\nConstruyendo el grafo con NetworkX...")
    # Forzamos que las columnas de los nodos sean de tipo 'string' ANTES de crear el grafo.
    # Esto asegura que pandas no las interprete como números.
    network_df['protein1'] = network_df['protein1'].astype(str)
    network_df['protein2'] = network_df['protein2'].astype(str)

    G = nx.from_pandas_edgelist(network_df, 'protein1', 'protein2', edge_attr='weight')

    # Verificar que las semillas CONVERTIDAS estén en la red
    original_seed_count = len(entrez_seeds)
    seed_genes_in_network = [seed for seed in entrez_seeds if G.has_node(seed)]
    print(
        f"Genes semilla (convertidos a Entrez) encontrados en la red: {len(seed_genes_in_network)} de {original_seed_count}.")

    if not seed_genes_in_network:
        print("Error: Ninguno de los genes semilla convertidos se encuentra en la red. Abortando.")
        return

    # --- Ejecución del algoritmo ---
    propagation_results = run_rwr(G, seed_genes_in_network)

    # --- Guardado de resultados ---
    if propagation_results is not None:
        print(f"Guardando los resultados en '{args.output}'...")
        output_dir = os.path.dirname(args.output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        propagation_results.to_csv(args.output, sep='\t', index=False)
        print("¡Análisis de propagación completado con éxito!")


if __name__ == "__main__":
    main()
