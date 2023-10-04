On considère deux blocs de masses respectives m1 et m2 liés l'un à l'autre par un ressort de constante de raideur k2 . Le bloc de masse m1 est lié à un point d'ancrage fixe par l'intermédiaire d'un ressort de constante de raideur k1 et, à l'autre extrémité du système, le bloc de masse m2 est lié à un point d'ancrage fixe par l'intermédiaire d'un ressort de constante de raideur k3 . La masse des ressorts est négligeable et on suppose que l'amplitude de déplacement des deux blocs est toujours suffisamment faible pour que la loi de Hooke soit vérifiée. Finalement, tous les frottements sont considérés comme négligeables.

1. quelle est la Matrice d'incidence $C$
2. quelle est la transposé Matrice d'incidence $C$ noté $C^T$
3. quelle est la matrice diagonale des masses $M$
4. quelle est la Matrice diagonale des raideur $D_K$
5. quelle est la matrice Matrice représentative des raideurs $K$ tel que $K = C^T .( D_K . C)$
6. quelle est la matrice du système $A$ tel que $M.\ddot{X} = A.X$
7. Ecrire l'équation du mouvement
8. Quelles sont les valeurs propre du système?
9. Quels sont les positions d'équilibres stable et instable du système ?
10. Quels sont les positions d'équilibre du système ?
11. Quels sont les positions d'equilibre stable et instable du systeme?
12. `Calcul du facteur de facteur de qualité Q ≡ ω0/(2λ)`
13. `y-a-t-il raisonnance? `
14. `Quelle est la fonction de transfert qui caractérise la réponse en fréquence de l’oscillateur étudié?`
15. `Calculer le gain du filtre.`
16. `Tracer Diagramme de Bode`
17. `Démontrer supermodularité`
18. `What is the difference between equilibrium and Accuracy`

---


1. **Matrice d'incidence $C$**
$C = \begin{bmatrix}-1 & 1 & 0\\0 & -1 & 1 \\\end{bmatrix}$

```python
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
```

2. **Transposé de la matrice $C$ noté $C^T$**
$$ C^T = \begin{bmatrix}
-1 & 0 \\
1 & -1 \\
0 & 1 \\
\end{bmatrix}$$

``` python
def get_transposed_incidence_matrix(incidence_matrix):
    """
    Get the transposed incidence matrix.
    Parameters:
    - incidence_matrix: The incidence matrix of the system.
    Returns the transposed incidence matrix(C^T).
    """
    return incidence_matrix.T
```

3. **Matrice des masses $M$**
$$
M = 
\begin{bmatrix}
m_1 & 0\\
0 & m_2\\
\end{bmatrix}
$$

```python
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
```
4. **Matrice diagonale des raideur $D_K$**

$$
D_K = 
\begin{bmatrix}
k_1 & 0 & 0 \\
0 & k_2 & 0 \\
0 & 0 & k_3 \\ 
\end{bmatrix}
$$

```python
def get_diagonal_stiffness_matrix(G):
    """
    Get the diagonal stiffness matrix of the system based on edge weights from a given graph.
    Parameters:
    - G: The input graph representing the physical system.
    Returns the diagonal stiffness matrix (D_K) .
    """
    raideurs = nx.get_edge_attributes(G, "k")
    return np.diag(list(raideurs.values()))
```

5. **Matrice representative des raideurs $K$**
on a : $K = C^T . D_K .C$

$$(1): D_K .C = 
\begin{bmatrix}
k_1 & 0 & 0 \\
0 & k_2 & 0 \\
0 & 0 & k_3 \\ 
\end{bmatrix}
.
\begin{bmatrix}
-1 & 0 \\
1 & -1 \\
0 & 1 \\
\end{bmatrix}
=
\begin{bmatrix}
-k_1 & k1-k2 & -k2 \\
0 & -k_2 & k2-k3 \\
0 & 0 & k_3 \\ 
\end{bmatrix}
$$
$$(2): C^T.(1) = 
\begin{bmatrix}
-1 & 1 & 0 \\
0 & -1 & 1 \\
\end{bmatrix}
.
\begin{bmatrix}
-k_1 & k1-k2 & -k2 \\
0 & -k_2 & k2-k3 \\
0 & 0 & k_3 \\ 
\end{bmatrix}
=
\begin{bmatrix}
k_1+k2 & -k2 \\
-k_2 & k2+k3 \\
\end{bmatrix}
$$
d'où:
$$K = C^T .( D_K .C) = 
\begin{bmatrix}
k_1+k_2 & -k_2 \\
-k_2 & k_2+k_3 \\
\end{bmatrix}$$
```python
def calculate_stiffness_matrix(incidence_matrix, diagonal_stiffness_matrix):
    """
    Calculate the stiffness matrix of the system.
    Parameters:
    - incidence_matrix: The incidence matrix of the system.
    - diagonal_stiffness_matrix: The diagonal stiffness matrix.
    Returns the stiffness matrix(K).
    """
    return np.dot(np.dot(incidence_matrix.T, diagonal_stiffness_matrix), incidence_matrix)
```


6. **Matrice representative du systeme $A$**
on a :     $A = M^{-1} \cdot K$ tel que $\ddot{X} = A \cdot X$
or: $M^{-1}$ inversible car M est une matrice diagonale.

On obtient alors:

$$A = 
\begin{bmatrix}
\frac{1}{m_1} & 0  \\
0 & \frac{1}{m_2}\\
\end{bmatrix}
.
\begin{bmatrix}
k_1+k_2 & -k_2  \\
-k_2 & k_2+k_3\\
\end{bmatrix}
$$
d'où:
$$ A = \begin{bmatrix}
\frac{k_1+k_2}{m1} & \frac{-k_2}{m1}  \\
\frac{-k_2}{m2} & \frac{k_2+k_3}{m2}\\
\end{bmatrix} $$
```python
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
```
1. On obtient l'équation du mouvement:
$$\ddot{X} = A.X$$
soit:
$$
\begin{bmatrix}
\dot{x}_1 \\
\dot{x}_2 \\
\end{bmatrix}
=
\begin{bmatrix}
\frac{k_1+k_2}{m1} & \frac{-k_2}{m1}  \\
\frac{-k_2}{m2} & \frac{k_2+k_3}{m2}\\
\end{bmatrix}
.
\begin{bmatrix}
x_1 \\
x_2 \\
\end{bmatrix}
$$

*Note: $A$ semble etre la matrice d'etat de notre systeme voir [[Retrouver Equation du neurone à partir des oscillateurs]]*

8. Quelles sont les valeurs propre du système?

Pour déterminer les valeurs propres λ de la matrice $A$, on doit résoudre l'équation caractéristique :$$det(A-\lambda I) = 0$$où $A$ est la matrice dont on veut trouver les valeurs propres, $λ$ est la variable propre, et $I$ est la matrice identité de la même dimension que $A$.
Dans notre cas, la matrice A est donnée par :
$$A = \begin{bmatrix}\frac{k_1+k_2}{m1} & \frac{-k_2}{m1}  \\\frac{-k_2}{m2} & \frac{k_2+k_3}{m2}\\\end{bmatrix}$$ 
Nous allons résoudre l'équation caractéristique pour cette matrice. Tout d'abord, définissons $I$ comme la matrice identité 2x2 :
$$I = \begin{bmatrix}1 & 0 \\0& 1\\\end{bmatrix}$$ 
Maintenant, formons la matrice $A−λI$:
$$A−λI = 
\begin{bmatrix}
\frac{k_1+k_2}{m1} & \frac{-k_2}{m1}  \\
\frac{-k_2}{m2} & \frac{k_2+k_3}{m2}\\
\end{bmatrix}
-
λ
\begin{bmatrix}1 & 0 \\0& 1\\\end{bmatrix}
=
\begin{bmatrix}
\frac{k_1+k_2}{m1}-λ  & \frac{-k_2}{m1}  \\
\frac{-k_2}{m2} & \frac{k_2+k_3}{m2}- λ  \\
\end{bmatrix}
$$

Maintenant, calculons le déterminant de cette matrice et égalons-le à zéro pour résoudre l'équation caractéristique :

$$det(A−λI) = 
\begin{vmatrix}
\begin{bmatrix}
\frac{k_1+k_2}{m1}-λ  & \frac{-k_2}{m1}  \\
\frac{-k_2}{m2} & \frac{k_2+k_3}{m2}- λ  \\
\end{bmatrix}
\end{vmatrix}$$


- Solve the systeme
- Implement Jacobian Solution

## Avec Frottement

**Matrice de frottement $B$**
Supposons que vous avez un vecteur des vitesses généralisées $\dot{X}$ (qui est le taux de changement des déplacements généralisés $X$). Les forces de frottement linéaires peuvent être exprimées sous forme matricielle comme suit :
$$F_{f1} = - B \dot{X}$$
où:
- $F_{f1}$ est le vecteur des forces de frottement.
- $B$ est une matrice diagonale de coefficients de frottement linéaire correspondant aux éléments de vitesse $\dot{x}$​.
Si on a $n$ degrés de liberté dans votre système, alors $B$ sera une matrice $n×n$ diagonale avec des coefficients de frottement sur la diagonale. Chaque coefficient de frottement correspond à la force de frottement associée à la vitesse du degré de liberté correspondant.
$$
B = 
\begin{vmatrix}
b_1 & 0 & 0 & \cdots & 0 \\
0 & b_2 & 0 & \cdots & 0 \\
0 & 0 & b_3 & \cdots & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & b_n
\end{vmatrix}
$$
```python
def calculate_damping_matrix(velocity_coefficient, mass_matrix):
    """
    Calculate the damping matrix of the system based on a damping coefficient and mass matrix.
    Parameters:
    - velocity_coefficient: The damping coefficient.
    - mass_matrix: The mass matrix of the system.
    Returns the damping matrix.
    """
    return velocity_coefficient * mass_matrix
```

**Equation du mouvement avec frottement**

$$\ddot{X} = A . X - B . \dot {X}$$
