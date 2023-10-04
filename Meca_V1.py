import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation
import sympy as sp


#==========================SHOW GRAPH
def show_graph(G):
    # Get info on nodes 
    node_name   = nx.get_node_attributes(G, "node_name")  # Obtient les noms des nœuds
    node_value  = nx.get_node_attributes(G, "node_value")  # get nodes values

    # Dessine le graphe en utilisant les valeurs des nœuds
    pos = nx.spring_layout(G, seed=42)  # Disposition des nœuds

    # Crée un dictionnaire d'étiquettes de nœuds en incluant les noms et les valeurs
    node_labels = {node: f"{node_name[node]}\n({node})" for node in G.nodes}

    #{node_name[node]}
    nx.draw(G, pos, labels=node_labels, node_size=500, node_color="skyblue", font_size=10, font_color="black", arrows=True)

    # Affiche les valeurs des arêtes à côté des arêtes
    edge_labels = {(u, v): poids for (u, v), poids in nx.get_edge_attributes(G, "k").items()}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10)

    plt.title("Graphe Orienté avec Valeurs sous les Noms des Nœuds")
    plt.show()

#==========================PERFORM QUICK ALGEBRAIC ANALYSIS
import numpy as np
def analyze_matrix(matrix):
    # Convert the matrix to a NumPy array
    arr = np.array(matrix)

    # Calculate the transpose of the matrix
    transpose = arr.T

    # Calculate the determinant of the matrix
    determinant = np.linalg.det(arr)

    # Calculate the inverse of the matrix
    try:
        inverse = np.linalg.inv(arr)
    except np.linalg.LinAlgError:
        inverse = None

    # Calculate the eigenvalues and eigenvectors of the matrix
    eigenvalues, eigenvectors = np.linalg.eig(arr)

    # Create a dictionary to store the results
    results = {
        'transpose': transpose.tolist(),
        'determinant': determinant,
        'inverse': inverse.tolist() if inverse is not None else None,
        'eigenvalues': eigenvalues.tolist(),
        'eigenvectors': eigenvectors.tolist()
    }
    print(f"Transpose:\n{results['transpose']}")
    print(f"Determinant: {results['determinant']}")
    print(f"Inverse:\n{results['inverse']}")
    print(f"Eigenvalues: {results['eigenvalues']}")
    print(f"Eigenvectors:\n{results['eigenvectors']}")
    return results


#==========================PERFORM QUICK GRAPH ANALYSIS

import networkx as nx
def analyze_graph(graph):
    # Get the number of nodes in the graph
    num_nodes = graph.number_of_nodes()

    # Get the number of edges in the graph
    num_edges = graph.number_of_edges()

    # Get the degree centrality of each node
    degree_centrality = nx.degree_centrality(graph)

    # Get the betweenness centrality of each node
    betweenness_centrality = nx.betweenness_centrality(graph)

    # Get the eigenvector centrality of each node
    eigenvector_centrality = nx.eigenvector_centrality(graph)

    # Get the clustering coefficient of each node
    clustering_coefficient = nx.clustering(graph)

    # Create a dictionary to store the results
    results = {
        'num_nodes': num_nodes,
        'num_edges': num_edges,
        'degree_centrality': degree_centrality,
        'betweenness_centrality': betweenness_centrality,
        'eigenvector_centrality': eigenvector_centrality,
        'clustering_coefficient': clustering_coefficient
    }

    print(f"Number of nodes: {results['num_nodes']}")
    print(f"Number of edges: {results['num_edges']}")
    print(f"Degree centrality: {results['degree_centrality']}")
    print(f"Betweenness centrality: {results['betweenness_centrality']}")
    print(f"Eigenvector centrality: {results['eigenvector_centrality']}")
    print(f"Clustering coefficient: {results['clustering_coefficient']}")
    return results

#==========================PERFORM MATRIX TO GRAPH
def matrix_to_graph(matrix):
    # Create an empty graph
    graph = nx.Graph()

    # Get the number of rows and columns in the matrix
    rows = len(matrix)
    cols = len(matrix[0])

    # Add nodes to the graph
    for i in range(rows):
        for j in range(cols):
            node = f"{i}{j}"
            graph.add_node(node)

    # Add edges to the graph
    for i in range(rows):
        for j in range(cols):
            node1 = f"{i}{j}"
            if j < cols - 1:
                node2 = f"{i}{j+1}"
                weight = matrix[i][j+1]
                graph.add_edge(node1, node2, weight=weight)
            if i < rows - 1:
                node3 = f"{i+1}{j}"
                weight = matrix[i+1][j]
                graph.add_edge(node1, node3, weight=weight)

    return graph

#========================== GET INCIDENCE MATRIX (ϕ)
def get_incidence_matrix(G):
    """
    Create an incidence matrix for the system based on the graph's edge connections.
    Parameters:
    - G: The input graph representing the physical system.
    Returns the incidence matrix (C) .
    """
    edges = list(G.edges)
    num_nodes = len(G.nodes)
    incidence_matrix = np.zeros((len(edges), num_nodes))
    
    for i, (source, target) in enumerate(edges):
        incidence_matrix[i, source] = -1
        incidence_matrix[i, target] = 1

    return incidence_matrix

#========================== GET TRANSPOSED MATRIX (ϕ_T)
def get_transposed_incidence_matrix(incidence_matrix):
    """
    Get the transposed incidence matrix.
    Parameters:
    - incidence_matrix: The incidence matrix of the system.
    Returns the transposed incidence matrix(C^T).
    """
    return incidence_matrix.T    

#========================== GET MASS MATRIX (M)
def get_mass_matrix(G):
    """
    Get the mass matrix of the system based on node values from a given graph.
    Parameters:
    - G: The input graph representing the physical system.
    Returns the diagonal mass matrix.
    """
    node_value = nx.get_node_attributes(G, "node_value")
    masses = list(node_value.values())
    return np.diag(masses)

#========================== GET MASS MATRIX (D_K)
def get_diagonal_stiffness_matrix(G):
    """
    Get the diagonal stiffness matrix of the system based on edge weights from a given graph.
    Parameters:
    - G: The input graph representing the physical system.
    Returns the diagonal stiffness matrix (D_K) .
    """
    raideurs = nx.get_edge_attributes(G, "k")
    return np.diag(list(raideurs.values()))    

#========================== GET STIFNESS MATRIX (K)
def calculate_stiffness_matrix(incidence_matrix, diagonal_stiffness_matrix):
    """
    Calculate the stiffness matrix of the system.
    Parameters:
    - incidence_matrix: The incidence matrix of the system.
    - diagonal_stiffness_matrix: The diagonal stiffness matrix.
    Returns the stiffness matrix(K).
    """
    return np.dot(np.dot(incidence_matrix.T, diagonal_stiffness_matrix), incidence_matrix)

#========================== GET STATE MATRIX (A)
def calculate_state_matrix(M, K, C):
    """
    Calculate the state matrix (A) of the system.
    Parameters:
    - M: The mass matrix of the system.
    - K: The stiffness matrix of the system.
    - C: The damping matrix of the system.
    Returns the state matrix (A).
    """
    if np.linalg.det(M) != 0:  # If the mass matrix is invertible
        M_inverse = np.linalg.inv(M)
        A = M_inverse @ K  # Add the damping matrix to the system
        return A
    else:
        print("\nMass Matrix M is not invertible. The system cannot be solved.")
        return None

#========================== GET OSCILLATION FREQUENCIES ✅
def calculate_oscillation_frequencies(M, K):
    """
    Calculate the oscillation frequencies of the system.

    Parameters:
    - M: The mass matrix of the system.
    - K: The stiffness matrix of the system.

    Returns a list of oscillation frequencies (in Hz).
    """
    # Solve the generalized eigenvalue problem to obtain eigenvalues (squared frequencies)
    eigenvalues = np.linalg.eigvals(np.linalg.inv(M).dot(K))

    # Calculate the oscillation frequencies by taking the square root of the eigenvalues
    frequencies = np.sqrt(np.abs(eigenvalues))

    # Sort the frequencies in ascending order
    frequencies.sort()

    # Convert frequencies from radians per second to Hertz
    frequencies_hz = frequencies / (2 * np.pi)

    return frequencies_hz


#========================== SOLVE SYSTEM DYNAMICS (G)
def solve_system_dynamics(G, initial_positions, initial_velocities, time_span, num_points):
    """
    Solve the dynamic behavior of a physical system represented as a graph.
    
    Parameters:
    - G: The input graph representing the physical system.
    - initial_positions: Initial positions of nodes.
    - initial_velocities: Initial velocities of nodes.
    - time_span: Time span for the simulation (e.g., [0, 10] for 0 to 10 seconds).
    - num_points: Number of time points in the simulation.
    
    Returns:
    - t: Time points at which the solution is computed.
    - X: Solution trajectories (positions) over time.
    """
    # Define a function representing the system's equations of motion
    def system_equations(t, X):
        # X contains positions and velocities as [x1, y1, x2, y2, ..., vx1, vy1, vx2, vy2, ...]
        positions = X[:len(X) // 2]
        velocities = X[len(X) // 2:]
        
        # Calculate the accelerations (second derivatives) from the matrix system
        accelerations = np.dot(A, positions)
        
        # Return the derivatives [dx1/dt, dy1/dt, dx2/dt, dy2/dt, ..., dvx1/dt, dvy1/dt, dvx2/dt, dvy2/dt, ...]
        return np.concatenate((velocities, accelerations))
    
    # Calculate the system's state matrix A based on the graph
    ϕ = get_incidence_matrix(G)
    D_K = get_diagonal_stiffness_matrix(G)
    K = calculate_stiffness_matrix(ϕ, D_K)
    M = get_mass_matrix(G)
    A = calculate_state_matrix(M, K, ϕ)
    
    # Combine initial positions and velocities into a single array
    initial_state = np.concatenate((initial_positions, initial_velocities))
    
    # Solve the system's dynamics using solve_ivp
    solution = solve_ivp(system_equations, time_span, initial_state, t_eval=np.linspace(*time_span, num_points))
    
    # Extract time points and solution trajectories
    t = solution.t
    X = solution.y
    
    return t, X 


#========================== FIND EQUILIBRIUM POSITION (G)
def find_equilibrium_position(G):
    """
    Find the equilibrium position of a physical system represented as a graph.
    
    Parameters:
    - G: The input graph representing the physical system.
    
    Returns:
    - equilibrium_positions: An array of equilibrium positions for each node.
    """
    # Calculate the system's state matrix A based on the graph
    ϕ = get_incidence_matrix(G)
    D_K = get_diagonal_stiffness_matrix(G)
    K = calculate_stiffness_matrix(ϕ, D_K)
    M = get_mass_matrix(G)
    A = calculate_state_matrix(M, K, ϕ)
    
    # Solve for the equilibrium positions by setting accelerations to zero
    num_nodes = len(G.nodes)
    zero_accelerations = np.zeros(num_nodes)
    equilibrium_positions = np.linalg.solve(A, zero_accelerations)
    
    return equilibrium_positions

#========================== FIND EQUILIBRIUM POSITION (A)
# Function to find the equilibrium positions
def find_equilibrium(A, gravity=0):
    """
    Find the equilibrium positions of a physical system represented as a graph.
    
    Parameters:
    - A: The state matrix of the system.
    - gravity: Gravitational acceleration (optional).
    
    Returns:
    - equilibrium_positions: An array of equilibrium positions for each node.
    """
    n = len(A)  # Number of nodes
    b = np.zeros(n)
    x = np.zeros(n)  # Define x to be an array of zeros
    for i in range(n):
        for j in range(n):
            b[i] -= A[i, j] * (0 - x[j])
        b[i] -= gravity
    x = np.linalg.solve(A, b)
    return x


#========================== FIND JACOBIAN EQUILIBRIUM POSITION (J)
def find_equilibrium_positions_Jacob(M, A):
    """
    Find the stable and unstable equilibrium positions of the system based on mass and coefficient matrices.
    Parameters:
    - M: The mass matrix of the system.
    - A: The coefficient matrix representing the equations of motion.
    Returns lists of stable and unstable equilibrium positions (Jacobian matrix J = M^-1 * A).
    """
    J = np.linalg.inv(M) @ A

    # Calculate the eigenvalues of the Jacobian matrix
    eigenvalues = np.linalg.eigvals(J)

    stable_equilibria = []
    unstable_equilibria = []

    for eigenvalue in eigenvalues:
        if np.real(eigenvalue) < 0:
            stable_equilibria.append(np.real(eigenvalue))
        elif np.real(eigenvalue) > 0:
            unstable_equilibria.append(np.real(eigenvalue))

    return stable_equilibria, unstable_equilibria


#========================== CREATION OF GRAPH (G)
def create_graph():
    """
    Create a new directed graph representing a physical system with nodes as masses and edges as springs.
    Returns the created graph.
    """
    G = nx.DiGraph()
    m1 = 100
    m2 = 2
    m3 = 1
    m4 = 10
    m5 = 5

    # Add 3 nodes to the graph with assigned values
    G.add_node(0, node_name="M0", node_value=m1)
    G.add_node(1, node_name="M1", node_value=m2)
    G.add_node(2, node_name="M2", node_value=m3)
    G.add_node(3, node_name="M3", node_value=m4)
    G.add_node(4, node_name="M4", node_value=m5)

    # Add weighted edges to the graph based on the stiffness matrix
    G.add_edge(0, 1, k = 1.0)
    G.add_edge(1, 0, k = 1.0)
    G.add_edge(1, 2, k = 2.0)
    G.add_edge(2, 1, k = 2.0)
    G.add_edge(3, 1, k = 3.0)
    G.add_edge(4, 1, k = 3.0)

    return G








def main():

    #----------------- TO DEFINE ----------------------------------------------
    initial_positions   = np.array([1.0, 2.0, 3.0, 4.0, 5.0])  # Initial positions
    initial_velocities  = np.array([1.0, 0.0, 0.0, 0.0, 0.0])  # Initial velocities
    time_span           = [0, 10]  # Time span for simulation
    num_points          = 100  # Number of time points
    #--------------------------------------------------------------------------

    G = create_graph()  # Assuming you have a function to create the graph

    ϕ = get_incidence_matrix(G)
    D_K = get_diagonal_stiffness_matrix(G)
    K = calculate_stiffness_matrix(ϕ, D_K)
    M = get_mass_matrix(G)
    A = calculate_state_matrix(M, K, ϕ)

    t, X = solve_system_dynamics(G, initial_positions, initial_velocities, time_span, num_points)


    equilibrium_positions_Test_1 = find_equilibrium_position(G)
    equilibrium_positions_Test_2 = find_equilibrium(A, gravity=0)
    stable, unstable = find_equilibrium_positions_Jacob(M, A)
    frequencies = calculate_oscillation_frequencies(M, K)
    

    show_graph(G)
    print("\nIncidence Matrix ϕ :\n", ϕ )
    print("\nDiagonal Stiffness Matrik D_K:\n",D_K)
    print("\nStiffness Matrix K:\n",K)
    print("\nMass MAtrix M:\n", M)
    print("\nState Matrix Representation A:\n",A)
    print("\nEquilibrium Positions:\n", equilibrium_positions_Test_1)
    print("\nEquilibrium Positions 2 option:\n", equilibrium_positions_Test_2)


    # ======== Stability Equilibrium position analysis 

    print("\nStable Equilibrium Positions:")
    for position in stable:
      print(position)

    print("\nUnstable Equilibrium Positions:")
    for position in range(len(unstable)):
        print("Node {}: {:.15f}".format(position+1, unstable[position]))

    # Step 9: Calculate oscillation frequencies

        # Print the frequencies
    print("\nOscillation Frequencies:")
    for i in range(len(frequencies)):
        print("Edge {}: {:.2f}".format(i+1, frequencies[i]))

if __name__ == "__main__":
    main()
