import networkx as nx
import numpy as np
from matplotlib import pyplot as plt

def draw_network(markov_matrix):
    num_states = len(markov_matrix)
    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes
    for state in range(num_states):
        G.add_node(state)

    # Add edges with weights
    for i in range(num_states):
        for j in range(i+1, num_states):  # Only consider upper triangular matrix
            if markov_matrix[i, j] > 0:
                G.add_edge(i, j, weight=markov_matrix[i, j], rad=0.2)  # Curved edge for upper triangular
            if markov_matrix[j, i] > 0:
                G.add_edge(j, i, weight=markov_matrix[j, i], rad=0.2) # Curved edge for lower triangular

    # Draw the graph
    pos = nx.circular_layout(G)


    lbls_from = {(i, j): f'{G[i][j]["weight"]*100:.2f}%' for i, j in G.edges if i<j}
    lbsl_to = {(i, j): f'{G[i][j]["weight"]*100:.2f}%' for i, j in G.edges if i>j}
    lbls = {**lbls_from, **lbsl_to}


    fig, ax = plt.subplots(figsize=(12, 9))

    params = {(0,1): (0.25, 0.45, '#FF6347',), (0,2): (0.25, 0.6, '#FFD700'), (2,1): (0.25, 0.6, '#20B2AA')}
    add_labels(G, lbls, params, ax, pos)

    colors = [params[k][2] for k in G.edges]

    G = nx.relabel_nodes(G, {0: 'Authoritarian', 1: 'Middle', 2: 'Democratic'})
    pos = nx.circular_layout(G)
    ax.set_title('Markov Chain of movements between different regiems')
    nx.draw(G, pos, with_labels=True,node_size=7000, node_color='skyblue', font_size=10, font_weight='bold', 
            width=[G[i][j]['weight']*150 for i, j in G.edges],connectionstyle='arc3,rad=0.4', arrowsize=50, 
            edge_color = colors, ax = ax)

    fig.show()
    


   
def add_labels(G, lbls, params,ax, pos):
    for k in list(params):
        params[(k[1], k[0])] = params[k]
    for edge in G.edges():
        a, loc, col = params[edge]
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        det = ((x1-x0)**2 + (y1-y0)**2)**0.5
        x = x0 + (x1-x0)*loc - a*(y0-y1)/det
        y = y0 + (y1 -y0)*loc - a*(x1-x0)/det
        ax.text(x-0.12, y, lbls[edge], fontsize=12,
            bbox=dict(facecolor='none', edgecolor='none'), 
            horizontalalignment='left')
        