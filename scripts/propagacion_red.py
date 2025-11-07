import pandas as pd
import networkx as nx
import numpy as np
import argparse
import os


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
    # 1. Preparar la matriz de transición normalizada por columnas
    nodes = list(graph.nodes())
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    adj_matrix = nx.to_numpy_array(graph, nodelist=nodes, weight='weight', dtype=float)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        degree = adj_matrix.sum(axis=0)
        # La matriz de transición M_ij es la probabilidad de pasar de j a i
        transition_matrix = adj_matrix / degree
        transition_matrix[np.isnan(transition_matrix)] = 0

    # 2. Inicializar los vectores de puntuación/probabilidad
    num_nodes = len(nodes)
    p0 = np.zeros(num_nodes) # Vector de reinicio
    for seed in seeds:
        if seed in node_to_idx:
            p0[node_to_idx[seed]] = 1.0
    p0 /= p0.sum() # Normalizar para que sume 1

    pt = np.copy(p0) # Vector de probabilidad en el tiempo t (iteración actual)

    # 3. Bucle iterativo del algoritmo RWR
    print("Iniciando Random Walk with Restart...")
    for i in range(max_iter):
        # Ecuación del RWR: p_{t+1} = (1 - alpha) * M * p_t + alpha * p_0
        pt_next = (1 - alpha) * np.dot(transition_matrix, pt) + alpha * p0
        
        # Comprobar la convergencia midiendo la diferencia entre vectores
        diff = np.linalg.norm(pt_next - pt, 1)
        print(f"Iteración {i+1}, Cambio: {diff:.6f}")
        
        if diff < tolerance:
            print(f"Convergencia alcanzada en la iteración {i+1}.")
            pt = pt_next
            break
        
        pt = pt_next

    # 4. Formatear y devolver los resultados
    results_df = pd.DataFrame({'gene': nodes, 'score': pt})
    return results_df.sort_values(by='score', ascending=False).reset_index(drop=True)

def main():
    """
    Función principal para parsear argumentos y ejecutar el flujo de trabajo.
    """
    parser = argparse.ArgumentParser(
        description="Realiza propagación en red con Random Walk with Restart (enfoque GUILD).",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--network", 
        required=True, 
        help="Ruta al archivo de la red. Formato esperado:\n"
             "  - 3 columnas (nodo1, nodo2, peso) separadas por tab/espacio.\n"
             "  - 2 columnas (nodo1, nodo2), se asumirá peso=1."
    )
    parser.add_argument(
        "--seeds", 
        required=True, 
        help="Ruta al archivo de genes semilla (un gen por línea sin cabecera)."
    )
    parser.add_argument(
        "--output", 
        required=True, 
        help="Ruta para guardar el archivo de resultados (formato TSV)."
    )
    parser.add_argument(
        "--alpha", 
        type=float, 
        default=0.85, 
        help="Probabilidad de reinicio (default: 0.85)."
    )
    
    args = parser.parse_args()

    # --- Carga de datos ---
    print(f"Cargando la red desde '{args.network}'...")
    try:
        # Intentar leer con 3 columnas, si falla, probar con 2
        try:
            network_df = pd.read_csv(args.network, sep='\s+', header=0)
        except ValueError:
            network_df = pd.read_csv(args.network, sep='\s+', header=None, names=['protein1_hugo', 'protein2_hugo'])
            network_df['weight'] = 1.0 # Asignar peso 1 si no se proporciona
            print("No se encontró columna de pesos, se asume un peso de 1.0 para todas las interacciones.")
    except FileNotFoundError:
        print(f"Error: No se pudo encontrar el archivo de red en '{args.network}'")
        return
    
    

    print(f"Cargando genes semilla desde '{args.seeds}'...")
    try:
        seeds_df = pd.read_csv(args.seeds, header=None, names=['gene'])
        seed_genes = seeds_df['gene'].tolist()
    except FileNotFoundError:
        print(f"Error: No se pudo encontrar el archivo de semillas en '{args.seeds}'")
        return

    # --- Construcción del Grafo ---
    print("Construyendo el grafo con NetworkX...")
    # Esta es la línea corregida
    G = nx.from_pandas_edgelist(network_df, 'protein1_hugo', 'protein2_hugo', edge_attr='combined_score')

    # Verificar que las semillas estén en la red
    original_seed_count = len(seed_genes)
    seed_genes_in_network = [seed for seed in seed_genes if G.has_node(seed)]
    print(f"Genes semilla encontrados en la red: {len(seed_genes_in_network)} de {original_seed_count}.")
    if not seed_genes_in_network:
        print("Error: Ninguno de los genes semilla se encuentra en la red. Abortando.")
        return

    # --- Ejecución del algoritmo ---
    propagation_results = run_rwr(G, seed_genes_in_network, alpha=args.alpha)

    # --- Guardado de resultados ---
    print(f"Guardando los resultados en '{args.output}'...")
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    propagation_results.to_csv(args.output, sep='\t', index=False)
    print("¡Análisis de propagación completado con éxito!")

if __name__ == "__main__":
    main()
