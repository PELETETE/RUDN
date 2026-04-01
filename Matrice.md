# Matrice de différentiation de Chebyshev
## 1. Principe général

La méthode de collocation de Chebyshev consiste à approcher une fonction \( y(x) \) par un polynôme d'interpolation de degré \( N \) aux points de Chebyshev. La **matrice de différentiation** \( D \) est alors l'opérateur linéaire discret qui associe au vecteur des valeurs de la fonction aux nœuds le vecteur des valeurs de sa dérivée aux mêmes nœuds.

Soient \( x_0, x_1, \dots, x_N \) les points de collocation (extrema de \( T_N \) ou racines de \( T_{N+1} \)). On note :

\[
\mathbf{y} = \begin{pmatrix} y(x_0) \\ y(x_1) \\ \vdots \\ y(x_N) \end{pmatrix}, \qquad
\mathbf{y}' = \begin{pmatrix} y'(x_0) \\ y'(x_1) \\ \vdots \\ y'(x_N) \end{pmatrix}
\]

La matrice de différentiation \( D \) d'ordre 1 satisfait :

\[
\mathbf{y}' = D \mathbf{y}
\]

De même, pour la dérivée seconde :

\[
\mathbf{y}'' = D^{(2)} \mathbf{y} = D^2 \mathbf{y}
\]

---

## 2. Choix des points de collocation

Deux choix classiques pour les points de Chebyshev :

### 2.1 Nœuds de Chebyshev-Gauss (racines de \( T_{N+1} \))

\[
x_k = \cos\left( \frac{2k+1}{2N+2} \pi \right), \quad k = 0, 1, \dots, N
\]

Ces points sont ouverts (n'incluent pas les bords \( \pm 1 \)).

### 2.2 Points de Chebyshev-Gauss-Lobatto (extrema de \( T_N \)) — **le plus utilisé en collocation**

\[
x_k = \cos\left( \frac{k\pi}{N} \right), \quad k = 0, 1, \dots, N
\]

Ces points incluent les bords \( x_0 = 1 \) et \( x_N = -1 \) (ou l'inverse selon l'ordre).  
C'est ce choix que nous adoptons car il simplifie l'imposition des conditions aux limites.

---

## 3. Construction de la matrice de différentiation

### 3.1 Méthode par interpolation polynomiale

À chaque point \( x_j \), on considère le polynôme d'interpolation de Lagrange de degré \( N \) passant par tous les points :

\[
y(x) \approx \sum_{j=0}^{N} y(x_j) \, \ell_j(x)
\]

où les polynômes de Lagrange sont :

\[
\ell_j(x) = \prod_{\substack{i=0 \\ i \ne j}}^{N} \frac{x - x_i}{x_j - x_i}
\]

La dérivée s'écrit :

\[
y'(x_k) = \sum_{j=0}^{N} y(x_j) \, \ell_j'(x_k)
\]

Ainsi, les coefficients de la matrice de différentiation sont :

\[
D_{kj} = \ell_j'(x_k)
\]

### 3.2 Formules explicites pour les points de Chebyshev-Gauss-Lobatto

Soient \( x_k = \cos(\theta_k) \) avec \( \theta_k = \frac{k\pi}{N} \), \( k = 0, \dots, N \).

On introduit les coefficients :

\[
c_k = \begin{cases}
2, & k = 0 \text{ ou } k = N \\
1, & 1 \le k \le N-1
\end{cases}
\]

Alors les éléments de la matrice \( D \) (pour la dérivée première) sont donnés par :

\[
D_{kj} = \frac{c_k}{c_j} \cdot \frac{(-1)^{k+j}}{x_k - x_j}, \quad k \ne j
\]

\[
D_{kk} = -\sum_{\substack{j=0 \\ j \ne k}}^{N} D_{kj}
\]

Plus explicitement, pour les termes diagonaux :

\[
D_{kk} = \begin{cases}
-\dfrac{2N^2 + 1}{6}, & k = 0 \\[8pt]
-\dfrac{x_k}{2(1 - x_k^2)}, & 1 \le k \le N-1 \\[8pt]
\dfrac{2N^2 + 1}{6}, & k = N
\end{cases}
\]

### 3.3 Matrice de dérivée seconde

La matrice de dérivée seconde peut être obtenue de deux façons :

1. **Par produit matriciel** : \( D^{(2)} = D \times D \) (attention : cela donne la dérivée seconde exacte dans l'espace polynomial, mais pas nécessairement la même qu'une matrice dédiée).

2. **Par formule directe** (plus stable numériquement pour les hauts degrés) :

\[
D^{(2)}_{kj} = \frac{2}{c_j} \cdot \frac{(-1)^{k+j}}{x_k - x_j} \left( \frac{x_k}{1-x_k^2} - \frac{x_j}{1-x_j^2} \right), \quad k \ne j
\]

\[
D^{(2)}_{kk} = -\sum_{\substack{j=0 \\ j \ne k}}^{N} D^{(2)}_{kj}
\]

Avec la même convention pour \( c_j \).

---

## 4. Algorithme de construction en pratique

Voici les étapes pour construire \( D \) et \( D^{(2)} \) :

```
1. Définir N (degré polynomial)
2. Créer les points x_k = cos(kπ/N) pour k = 0..N
3. Définir les coefficients c_k : c_0 = c_N = 2, c_k = 1 sinon
4. Pour chaque k = 0..N
     Pour chaque j = 0..N
         Si k ≠ j
             D_{kj} = (c_k/c_j) * ((-1)^{k+j}) / (x_k - x_j)
         Sinon
             D_{kk} = -sum_{j≠k} D_{kj}
5. Calculer D2 = D * D  (produit matriciel classique)
   OU utiliser les formules directes pour D2
```

---

## 5. Transformation sur un intervalle [a, b]

Si le problème est posé sur \( [a, b] \) au lieu de \([-1, 1]\), on utilise la transformation affine :

\[
\xi = \frac{b-a}{2} x + \frac{a+b}{2}, \quad x \in [-1, 1]
\]

Les matrices de différentiation se transforment comme :

\[
D_{[a,b]} = \frac{2}{b-a} \, D_{[-1,1]}
\]

\[
D^{(2)}_{[a,b]} = \left( \frac{2}{b-a} \right)^2 D^{(2)}_{[-1,1]}
\]

---

## 6. Exemple pour N = 4 (points de Chebyshev-Gauss-Lobatto)

Points : \( x_0 = 1, x_1 = \cos(\pi/4) = \sqrt{2}/2 \approx 0.7071, x_2 = 0, x_3 = -\sqrt{2}/2 \approx -0.7071, x_4 = -1 \)

Matrice \( D \) (dérivée première) pour \( N = 4 \) :

\[
D \approx \begin{pmatrix}
-6.6667 & 8.0000 & -2.6667 & 1.3333 & -0.6667 \\
-1.4142 & 0.0000 & 1.4142 & -0.7071 & 0.7071 \\
0.6667 & -2.6667 & 0.0000 & 2.6667 & -0.6667 \\
-0.7071 & 0.7071 & -1.4142 & 0.0000 & 1.4142 \\
0.6667 & -1.3333 & 2.6667 & -8.0000 & 6.6667
\end{pmatrix}
\]

(Noter l'antisymétrie quasi-évidente et les valeurs élevées aux bords.)

---

## 7. Propriétés importantes

| Propriété | Description |
|-----------|-------------|
| **Exactitude polynomiale** | Pour tout polynôme de degré ≤ N, \( D \mathbf{y} \) donne exactement la dérivée aux nœuds. |
| **Matrice pleine** | Contrairement aux différences finies, \( D \) est dense → coût \( O(N^2) \) pour l'application, \( O(N^3) \) pour l'inversion. |
| **Stabilité** | Bien conditionnée pour des modérés, mais le conditionnement croît comme \( O(N^2) \). |
| **Convergence spectrale** | Pour les fonctions lisses, l'erreur de dérivation décroît exponentiellement avec \( N \). |

---

## 8. Utilisation pour une ODE du second ordre

Pour une équation du type :

\[
y''(x) + p(x) y'(x) + q(x) y(x) = f(x)
\]

Avec conditions aux limites, par exemple :

\[
y(-1) = \alpha, \quad y(1) = \beta
\]

Le système discrétisé s'écrit :

\[
\left( D^{(2)} + \operatorname{diag}(p(x_i)) \cdot D + \operatorname{diag}(q(x_i)) \right) \mathbf{y} = \mathbf{f}
\]

Les équations correspondant aux points de bord sont remplacées par les conditions aux limites :

\[
y_0 = \alpha, \quad y_N = \beta
\]

On obtient ainsi un système linéaire \( (N+1) \times (N+1) \) à résoudre.
