
# Tarea 2: Network Propagation

La **propagación en redes** es una técnica computacional versátil, ampliamente utilizada en ciencias biológicas para analizar estructuras de redes. Su idea fundamental consiste en difundir "señales" a través de una red, donde los nodos representan elementos como genes o proteínas, y las aristas sus interacciones o correlaciones. Este proceso amplifica señales débiles partiendo de un subconjunto inicial de nodos conocidos, revelando patrones y relaciones no evidentes en interacciones directas.

Una característica clave es su capacidad para considerar estructuras locales y globales: la **información se propaga iterativamente por nodos conectados**, identificando elementos relevantes en vías más amplias. Es útil para **integrar datos multi-ómicos**, **predecir funciones proteicas** y **descubrir relaciones en redes biológicas** a gran escala.

---

## Objetivo
La idea general del proyecto es que al tomar una lista de genes iniciales (llamados "semillas"), utilizando una red de interacciones entre proteínas/genes, se encuentren otros genes que estén funcionalmente relacionados o "cerca" de las semillas en la red.

A través de un un script en Python se va a implementar un ejemplo de propagación en redes utilizando algún algoritmo de **GUILD** y/o **DIAMOnD**. 

La semilla empleada será el conjunto de genes **ENO1**, **PGK1** y **HK2**. Se proporcionan archivos de red formateada para ambos algoritmos, así como una red filtrada obtenida de STRING y un script de ejemplo para procesarla.

Este script implementa un **análisis de propagación en red** utilizando el **algoritmo Random Walk with Restart (RWR)**. Este es el mismo enfoque algorítmico que utilizan herramientas bioinformáticas de referencia como **GUILD** para **priorizar genes**. Al implementar RWR, nuestro script reproduce la funcionalidad central del método de puntuación de GUILD.

## Herramienta GUILD
GUILD (***Genes Underlying Inheritance Disorders***) es un paquete de software bioinformático publicado en un artículo científico. Su propósito es tomar una lista de genes semilla y una red de interacciones para **priorizar** otros **genes** que podrían estar implicados en el mismo proceso biológico o enfermedad.

---

## Implementación del script

El script de propagación (`propagacion_red.py`) no incluye una etapa de conversión de identificadores a HUGO por diseño. El script se enfoca exclusivamente en la fase de análisis algorítmico, asumiendo como prerrequisito que los datos de entrada ya han sido pre-procesados y estandarizados. Esta premisa se cumple, ya que tanto los archivos de red proporcionados (`string_network_filtered_hugo-400.tsv`) como el archivo de semillas (`genes_seed.txt`) utilizan consistentemente la nomenclatura de símbolos HUGO. De esta forma se separan las etapas de pre-procesamiento de datos y análisis en cuestión, con el fin de generar eficiencia en el código.

### Algoritmo RWR: `run_rwr`

***Random Walk with Restart*** es un algoritmo matemático de propagación en redes que simula un caminante aleatorio en un grafo, donde los nodos representan entidades como genes o proteínas y las aristas sus interacciones. Se basa en calcular la "proximidad" de nodos en cualquier grafo y tiene como exclusividad que incluye un **mecanismo de reinicio** que equilibra influencias locales y globales.

1. **Inicio**: Se asigna una puntuación inicial a los genes semilla.
2. **Propagación**: En cada paso, esta puntuación se "propaga" desde cada gen a sus vecinos en la red.
3. **Reinicio**: En cada paso, también hay una probabilidad de que la puntuación "salte" de vuelta a los nodos semilla originales. Esto asegura que los genes más cercanos a las semillas reciban una mayor puntuación.
4. **Convergencia**: El proceso se repite hasta que las puntuaciones de todos los genes en la red se estabilizan.

## Resultado

El resultado final será una **lista de todos los genes de la , clasificados por su puntuación final**, que valora la proximidad funcional a las semillas. Los genes con las puntuaciones más altas (después de las semillas) son los candidatos más prometedores para estar funcionalmente relacionados con el proceso biológico de los genes semilla.


### Características de la implementación

1. Será ejecutable desde la línea de comandos (CLI).
2. Aceptará como argumentos la ruta al archivo de la red, la ruta al archivo de los genes semilla y la ruta del archivo de salida.
3. Cargará la red y las semillas usando las librerías pandas y networkx.
4. Implementará el algoritmo RWR.
5. Guardará los resultados en un archivo de texto, ordenados por puntuación.

## Estructura del repositorio

```

/network\_propagation/
├── data/
│   ├── network\_guild.txt                        # Red formateada para GUILD
│   ├── network\_diamond.txt                      # Red formateada para DIAMOnD
│   ├── string\_network\_filtered\_hugo-400.tsv     # Red filtrada de STRING
│   └── genes\_seed.txt                           # Genes semilla: ENO1, PGK1, HK2
├── scripts/
│   ├── process\_STRING.py                        # Script de ejemplo para procesar la red
│   └── propagacion_red.py                        # Script creado
├── README.md                                    # Este archivo
└── requirements.txt                             # Dependencias: networkx, pandas

```

## Dependencias recomendadas

Incluye en `requirements.txt` las librerías necesarias para ejecutar tu script. Por ejemplo:

```
pandas
networkx
numpy
argparse
os
mygene
```

