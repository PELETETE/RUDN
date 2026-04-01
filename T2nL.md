La méthode de collocation de Tchebychev constitue ainsi un outil de choix pour la résolution numérique d'EDO, alliant puissance théorique et facilité d'implémentation pratique.






# Rapport de Recherche
## Résolution d'équations différentielles ordinaires du second ordre par méthodes de collocation de Tchebychev

---

## Introduction : Fondements des méthodes spectrales

Les méthodes spectrales constituent une classe puissante de techniques numériques pour la résolution d'équations différentielles. Contrairement aux méthodes locales (différences finies, éléments finis) qui n'utilisent qu'une information locale autour de chaque point, les méthodes spectrales représentent la solution comme une combinaison linéaire de fonctions globales, définies sur tout le domaine.

### Critères de choix d'une base spectrale

Le choix de la base de fonctions est crucial et repose sur plusieurs critères :

| Critère | Importance |
|---------|------------|
| **Convergence rapide** | Une bonne base doit permettre une convergence exponentielle (dite "spectrale") avec le nombre de termes |
| **Facilité de dérivation/intégration** | Les opérateurs différentiels doivent s'appliquer simplement aux fonctions de base |
| **Orthogonalité** | Une base orthogonale simplifie considérablement les calculs et améliore le conditionnement |
| **Complétude** | Toute fonction raisonnable doit pouvoir être approchée dans la base |

### Problèmes non-périodiques et choix des polynômes de Tchebychev

Pour les problèmes périodiques, les séries de Fourier constituent un choix naturel. Cependant, pour les problèmes non-périodiques sur un intervalle fini, l'utilisation de points équidistants avec des polynômes conduit au phénomène de Runge : des oscillations divergentes apparaissent aux bords de l'intervalle.

Les **polynômes de Tchebychev** offrent une solution élégante à ce problème. Leur distribution particulière des racines (plus denses aux bords) minimise l'erreur d'interpolation et évite le phénomène de Runge. Cette propriété, appelée **minimax**, garantit que l'erreur maximale est aussi petite que possible.

---

## Section 1 : Méthodes de calcul des coefficients spectraux

L'objectif fondamental des méthodes spectrales est de transformer une équation différentielle continue en un système algébrique discret. Pour une équation $Lu = f$ sur un domaine $\Omega$, on approxime $u$ par une somme finie :

$$u_N(x) = \sum_{k=0}^{N} \hat{u}_k \phi_k(x)$$

Le résidu est défini par $R(x; \hat{u}) = Lu_N - f$. Différentes stratégies existent pour minimiser ce résidu.

### 1.1. Méthode de Galerkin

La méthode de Galerkin est l'approche la plus "pure" mathématiquement. On force le résidu à être orthogonal à chaque fonction de la base :

$$\langle R, \phi_j \rangle = \int_{\Omega} R(x) \phi_j(x) w(x) dx = 0 \quad \text{pour } j = 0, \dots, N$$

**Condition essentielle** : Chaque fonction de base $\phi_k$ doit satisfaire individuellement les conditions aux limites.

**Avantages** :
- Conduit à des matrices symétriques et bien conditionnées
- Fondement mathématique solide

**Inconvénients** :
- Construction d'une base adaptée aux conditions aux limites parfois complexe
- Gestion très difficile des non-linéarités

### 1.2. Méthode Tau de Lanczos

La méthode Tau étend Galerkin aux cas où les fonctions de base ne satisfont pas les conditions aux limites. L'idée est d'ajouter un terme correcteur ("Tau") au résidu.

**Avantages** :
- Travail direct dans l'espace des coefficients
- Pas besoin de construire une base spéciale

**Inconvénients** :
- Implémentation plus délicate
- Gestion extrêmement complexe des non-linéarités (produits de convolution dans l'espace spectral)

### 1.3. Méthode de Collocation (Pseudospectrale)

La méthode de collocation est la plus intuitive et la plus utilisée en pratique, **particulièrement pour les problèmes non-linéaires**. On exige que le résidu soit strictement nul en des points précis, appelés **nœuds de collocation** :

$$R(x_j) = 0 \quad \text{pour } j = 0, \dots, N$$

Cela revient à tester le résidu avec des fonctions Delta de Dirac $\delta(x - x_j)$. Pour une EDO non-linéaire du second ordre :

$$y''(x) = f(x, y, y')$$

la discrétisation donne un système non-linéaire :

$$D^2 \mathbf{y} = \mathbf{f}(x, \mathbf{y}, D\mathbf{y})$$

**Avantages décisifs pour les non-linéarités** :
- Pas d'intégrales complexes à calculer
- Les non-linéarités sont évaluées simplement aux nœuds
- S'adapte naturellement à la méthode de Newton
- Implémentation directe dans l'espace physique

---

## Section 2 : Propriétés des polynômes de Tchebychev

# Polynômes de Chebyshev : définitions et propriétés
## 1. Polynômes de Chebyshev de première espèce \( T_n(x) \)

### 1.1 Définition trigonométrique

Pour \( n \in \mathbb{N} \) et \( x \in [-1, 1] \) :

\[
T_n(x) = \cos(n \arccos x)
\]

Avec le changement de variable \( \theta = \arccos x \), on a \( x = \cos \theta \) et :

\[
T_n(\cos \theta) = \cos(n\theta), \quad \theta \in [0, \pi]
\]

### 1.2 Définition par récurrence

\[
\begin{cases}
T_0(x) = 1 \\[4pt]
T_1(x) = x \\[4pt]
T_{n+1}(x) = 2x\,T_n(x) - T_{n-1}(x), \quad n \ge 1
\end{cases}
\]

### 1.3 Expression explicite

Pour \( x \in \mathbb{R} \) (extension analytique) :

\[
T_n(x) = \frac{1}{2} \left[ \left( x + \sqrt{x^2 - 1} \right)^n + \left( x - \sqrt{x^2 - 1} \right)^n \right]
\]

Ou sous forme polynomiale :

\[
T_n(x) = \frac{n}{2} \sum_{k=0}^{\lfloor n/2 \rfloor} (-1)^k \frac{(n-k-1)!}{k! \, (n-2k)!} (2x)^{n-2k}
\]

Plus explicitement, les premiers termes :

\[
\begin{aligned}
T_0(x) &= 1 \\
T_1(x) &= x \\
T_2(x) &= 2x^2 - 1 \\
T_3(x) &= 4x^3 - 3x \\
T_4(x) &= 8x^4 - 8x^2 + 1 \\
T_5(x) &= 16x^5 - 20x^3 + 5x
\end{aligned}
\]

### 1.4 Orthogonalité

Sur l'intervalle \([-1, 1]\) avec le poids \( w(x) = \dfrac{1}{\sqrt{1-x^2}} \) :

\[
\int_{-1}^{1} T_m(x) \, T_n(x) \, \frac{dx}{\sqrt{1-x^2}} =
\begin{cases}
0, & m \ne n \\[4pt]
\pi, & m = n = 0 \\[4pt]
\dfrac{\pi}{2}, & m = n \ge 1
\end{cases}
\]

### 1.5 Points de collocation

**Nœuds (racines)** : \( n \) racines de \( T_n(x) \) dans \((-1, 1)\) :

\[
x_k = \cos\left( \frac{2k-1}{2n} \pi \right), \quad k = 1, 2, \dots, n
\]

**Extrema** : \( n+1 \) points incluant les bords :

\[
x_k = \cos\left( \frac{k\pi}{n} \right), \quad k = 0, 1, \dots, n
\]

---

## 2. Polynômes de Chebyshev de deuxième espèce \( U_n(x) \)

### 2.1 Définition trigonométrique

Pour \( n \in \mathbb{N} \) et \( x \in [-1, 1] \) :

\[
U_n(x) = \frac{\sin\big( (n+1) \arccos x \big)}{\sin(\arccos x)} = \frac{\sin\big( (n+1)\theta \big)}{\sin \theta}
\]

avec \( x = \cos \theta \).

### 2.2 Définition par récurrence

\[
\begin{cases}
U_0(x) = 1 \\[4pt]
U_1(x) = 2x \\[4pt]
U_{n+1}(x) = 2x\,U_n(x) - U_{n-1}(x), \quad n \ge 1
\end{cases}
\]

### 2.3 Expression explicite

\[
U_n(x) = \sum_{k=0}^{\lfloor n/2 \rfloor} (-1)^k \binom{n-k}{k} (2x)^{n-2k}
\]

Premiers termes :

\[
\begin{aligned}
U_0(x) &= 1 \\
U_1(x) &= 2x \\
U_2(x) &= 4x^2 - 1 \\
U_3(x) &= 8x^3 - 4x \\
U_4(x) &= 16x^4 - 12x^2 + 1 \\
U_5(x) &= 32x^5 - 32x^3 + 6x
\end{aligned}
\]

### 2.4 Orthogonalité

Sur \([-1, 1]\) avec le poids \( w(x) = \sqrt{1-x^2} \) :

\[
\int_{-1}^{1} U_m(x) \, U_n(x) \, \sqrt{1-x^2} \, dx =
\begin{cases}
0, & m \ne n \\[4pt]
\dfrac{\pi}{2}, & m = n
\end{cases}
\]

---

## 3. Relations entre \( T_n \) et \( U_n \)

### 3.1 Dérivées

\[
\frac{d}{dx} T_n(x) = n \, U_{n-1}(x)
\]

\[
\frac{d}{dx} U_n(x) = \frac{(n+1) T_{n+1}(x) - x U_n(x)}{x^2 - 1} \quad \text{(pour } x \ne \pm 1\text{)}
\]

Plus simplement :

\[
\frac{d}{dx} U_n(x) = \frac{1}{1-x^2} \big( (n+1) T_{n+1}(x) - x U_n(x) \big)
\]

### 3.2 Relations algébriques

\[
T_n(x) = U_n(x) - x U_{n-1}(x)
\]

\[
(1-x^2) U_{n-1}(x) = x T_n(x) - T_{n+1}(x)
\]

\[
T_n'(x) = n U_{n-1}(x)
\]

---

## 4. Transformation sur un intervalle \([a, b]\)

Pour un problème sur \([a, b]\), on utilise la transformation affine :

\[
x = \frac{2\xi - (a+b)}{b-a} \quad \Longleftrightarrow \quad \xi = \frac{b-a}{2} x + \frac{a+b}{2}
\]

où \( x \in [-1, 1] \) est la variable de Chebyshev et \( \xi \in [a, b] \) la variable physique.

Les polynômes deviennent :

\[
T_n^{(a,b)}(\xi) = T_n\!\left( \frac{2\xi - (a+b)}{b-a} \right)
\]

Les points de collocation dans \([a, b]\) sont alors :

\[
\xi_k = \frac{b-a}{2} \cos\!\left( \frac{k\pi}{N} \right) + \frac{a+b}{2}, \quad k = 0, \dots, N
\]

## Section 3 : Opérateurs spectraux et résolution d'EDO non-linéaires

### 3 Operateurs spectraux: Explication bien dosee.

# 1. Principe général

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
Nous suposons ici donc qu'une construction basee sur de telles regles mathemathematiques doit en 
retour nous fournir un environnement de d'estimation suivant les proprietes resumees en ce petit tableau ci-dessous.

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




### 3.1. Matrice de différenciation de Tchebychev(Resume en vu de comprendre accompagnee des comparaisons)

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.linalg import solve

def cheb_matrix(N):
    """
    Calcule la matrice de différenciation de Tchebychev D
    et les nœuds de collocation x pour N points
    """
    if N == 0: 
        return 0, np.array([1.0])
    
    # Nœuds de Tchebychev-Gauss-Lobatto
    x = np.cos(np.pi * np.arange(N + 1) / N)
    
    # Coefficients de pondération
    c = np.ones(N + 1)
    c[0], c[-1] = 2, 2
    c = c * (-1)**np.arange(N + 1)
    
    # Construction de la matrice
    X = np.tile(x, (N + 1, 1))
    dX = X - X.T
    
    # Éviter la division par zéro sur la diagonale
    D = np.outer(c, 1/c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    
    return D, x
```

### 3.2. Méthode de Newton pour systèmes non-linéaires

Pour un système non-linéaire $F(\mathbf{y}) = 0$, la méthode de Newton itère :

$$\mathbf{y}^{(k+1)} = \mathbf{y}^{(k)} - J(\mathbf{y}^{(k)})^{-1} F(\mathbf{y}^{(k)})$$

où $J$ est la matrice jacobienne.

### 3.3. Exemple 1 : Équation de Bratu (non-linéarité exponentielle)

L'équation de Bratu modélise la combustion et présente une bifurcation :
$$y''(x) + \lambda e^{y(x)} = 0, \quad y(-1) = y(1) = 0$$

```python
def bratu_nonlinear(N, lambda_param):
    """
    Résout l'équation de Bratu: y'' + lambda * exp(y) = 0, y(±1)=0
    par collocation de Tchebychev et méthode de Newton
    """
    # Matrices de différenciation
    D, x = cheb_matrix(N)
    D2 = D @ D
    
    # Initialisation (solution nulle)
    y = np.zeros(N + 1)
    
    # Paramètres de Newton
    max_iter = 50
    tol = 1e-12
    
    for iteration in range(max_iter):
        # Calcul de la fonction résiduelle F(y)
        # F(y) = D2 @ y + lambda * exp(y)
        F = D2 @ y + lambda_param * np.exp(y)
        
        # Application des conditions aux limites
        # On préserve les équations aux bords pour Newton
        # mais on forcera y[0]=0 et y[-1]=0 après correction
        
        # Calcul de la matrice jacobienne J
        # J = D2 + lambda * diag(exp(y))
        J = D2.copy()
        for i in range(N + 1):
            J[i, i] += lambda_param * np.exp(y[i])
        
        # Imposition des conditions aux limites dans le jacobien
        # y(1) = 0 -> ligne 0
        J[0, :] = 0
        J[0, 0] = 1
        F[0] = y[0] - 0
        
        # y(-1) = 0 -> ligne N
        J[-1, :] = 0
        J[-1, -1] = 1
        F[-1] = y[-1] - 0
        
        # Résolution du système linéaire J * delta = -F
        delta = solve(J, -F)
        
        # Mise à jour de la solution
        y = y + delta
        
        # Vérification de la convergence
        if np.linalg.norm(delta, np.inf) < tol:
            print(f"Convergence atteinte en {iteration+1} itérations")
            break
    
    return x, y

# Test avec différentes valeurs de lambda
plt.figure(figsize=(12, 8))

for lambda_val in [1.0, 2.0, 3.0, 3.5]:
    x, y = bratu_nonlinear(30, lambda_val)
    plt.plot(x, y, 'o-', label=f'λ = {lambda_val}')

plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.title('Équation de Bratu: $y\'\' + \\lambda e^{y} = 0$', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

### 3.4. Exemple 2 : Équation de Duffing (oscillateur non-linéaire)

L'équation de Duffing modélise un oscillateur avec raideur non-linéaire :
$$y''(x) + \omega^2 y(x) + \varepsilon y^3(x) = 0, \quad y(-1) = 0, \quad y(1) = 1$$

```python
def duffing_nonlinear(N, omega, epsilon):
    """
    Résout l'équation de Duffing: y'' + omega^2 y + epsilon y^3 = 0
    avec y(-1)=0, y(1)=1
    """
    D, x = cheb_matrix(N)
    D2 = D @ D
    
    # Initialisation (solution linéaire)
    y = np.linspace(0, 1, N + 1)
    
    max_iter = 50
    tol = 1e-12
    
    for iteration in range(max_iter):
        # Fonction résiduelle
        F = D2 @ y + omega**2 * y + epsilon * y**3
        
        # Jacobienne
        J = D2.copy()
        for i in range(N + 1):
            J[i, i] += omega**2 + 3 * epsilon * y[i]**2
        
        # Conditions aux limites
        # y(1) = 1 (x[0] = 1)
        J[0, :] = 0
        J[0, 0] = 1
        F[0] = y[0] - 1
        
        # y(-1) = 0 (x[-1] = -1)
        J[-1, :] = 0
        J[-1, -1] = 1
        F[-1] = y[-1] - 0
        
        # Mise à jour
        delta = solve(J, -F)
        y = y + delta
        
        if np.linalg.norm(delta, np.inf) < tol:
            print(f"Convergence Duffing en {iteration+1} itérations")
            break
    
    return x, y

# Comparaison des effets non-linéaires
plt.figure(figsize=(12, 8))

x_linear, y_linear = duffing_nonlinear(30, omega=5.0, epsilon=0.0)
plt.plot(x_linear, y_linear, 'b-', linewidth=2, label='Linéaire (ε=0)')

x_nonlin, y_nonlin = duffing_nonlinear(30, omega=5.0, epsilon=10.0)
plt.plot(x_nonlin, y_nonlin, 'r-', linewidth=2, label='Non-linéaire (ε=10)')

plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.title('Équation de Duffing: $y\'\' + \\omega^2 y + \\varepsilon y^3 = 0$', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

### 3.5. Exemple 3 : Équation avec non-linéarité quadratique et conditions de Neumann

Considérons maintenant un problème mixte avec conditions de Neumann :
$$y''(x) + y'(x) + y^2(x) = \cos(\pi x), \quad y'(-1) = 0, \quad y(1) = 0$$

```python
def mixte_nonlinear(N):
    """
    Résout: y'' + y' + y^2 = cos(πx)
    avec y'(-1)=0, y(1)=0
    """
    D, x = cheb_matrix(N)
    D2 = D @ D
    
    # Initialisation
    y = np.zeros(N + 1)
    
    max_iter = 50
    tol = 1e-12
    
    # Second membre
    f = np.cos(np.pi * x)
    
    for iteration in range(max_iter):
        # Résidu: R = y'' + y' + y^2 - f
        R = D2 @ y + D @ y + y**2 - f
        
        # Jacobienne: J = D2 + D + 2*diag(y)
        J = D2 + D
        for i in range(N + 1):
            J[i, i] += 2 * y[i]
        
        # Condition de Neumann à x = -1 (dernière ligne)
        # y'(-1) = 0 → D[-1, :] @ y = 0
        J[-1, :] = D[-1, :]
        R[-1] = D[-1, :] @ y - 0
        
        # Condition de Dirichlet à x = 1 (première ligne)
        # y(1) = 0
        J[0, :] = 0
        J[0, 0] = 1
        R[0] = y[0] - 0
        
        # Mise à jour Newton
        delta = solve(J, -R)
        y = y + delta
        
        if np.linalg.norm(delta, np.inf) < tol:
            print(f"Convergence en {iteration+1} itérations")
            break
    
    return x, y

# Résolution et visualisation
x, y = mixte_nonlinear(25)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(x, y, 'ro-', markersize=6)
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.title('Solution $y(x)$', fontsize=14)
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
# Calcul de y' pour vérifier la condition de Neumann
D, _ = cheb_matrix(25)
yp = D @ y
plt.plot(x, yp, 'b-o', markersize=4)
plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
plt.xlabel('x', fontsize=12)
plt.ylabel("y'(x)", fontsize=12)
plt.title('Dérivée $y\'(x)$', fontsize=14)
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
```

### 3.6. Implémentation générique pour EDO non-linéaire

Voici un solveur générique capable de traiter toute EDO non-linéaire du second ordre :

```python
def solveur_nonlinear_generique(N, f_nonlinear, conditions_limites, 
                                 y_init=None, max_iter=100, tol=1e-12):
    """
    Solveur générique pour EDO non-linéaire du second ordre
    y'' = f(x, y, y')
    
    Paramètres:
    - N: nombre de points
    - f_nonlinear: fonction f(x, y, yp) retournant la non-linéarité
    - conditions_limites: liste de tuples (type, position, valeur)
      type: 'D' pour Dirichlet, 'N' pour Neumann
      position: 0 (x=1) ou -1 (x=-1)
    - y_init: condition initiale (None = zéro)
    """
    D, x = cheb_matrix(N)
    D2 = D @ D
    
    # Initialisation
    if y_init is None:
        y = np.zeros(N + 1)
    else:
        y = y_init.copy()
    
    for iteration in range(max_iter):
        # Calcul des dérivées
        yp = D @ y
        ypp = D2 @ y
        
        # Résidu: R = ypp - f(x, y, yp)
        R = ypp - f_nonlinear(x, y, yp)
        
        # Jacobienne: J = D2 - df/dy - df/dyp * D
        # Approximation par différences finies pour la généralité
        J = D2.copy()
        epsilon = 1e-8
        
        for i in range(N + 1):
            # Perturbation pour df/dy
            y_pert = y.copy()
            y_pert[i] += epsilon
            yp_pert = D @ y_pert
            f_pert = f_nonlinear(x[i], y_pert[i], yp_pert[i])
            
            y_current = y[i]
            yp_current = yp[i]
            f_current = f_nonlinear(x[i], y_current, yp_current)
            
            df_dy = (f_pert - f_current) / epsilon
            
            # Perturbation pour df/dyp
            yp_pert2 = yp.copy()
            yp_pert2[i] += epsilon
            # Il faudrait perturber y pour affecter yp de manière cohérente
            # Version simplifiée:
            df_dyp = 0  # À améliorer selon le problème
            
            J[i, i] += -df_dy  # Contribution de -df/dy
            
        # Application des conditions aux limites
        for cond_type, pos, val in conditions_limites:
            if cond_type == 'D':  # Dirichlet
                J[pos, :] = 0
                J[pos, pos] = 1
                R[pos] = y[pos] - val
            elif cond_type == 'N':  # Neumann
                J[pos, :] = D[pos, :]
                R[pos] = yp[pos] - val
        
        # Mise à jour Newton
        try:
            delta = solve(J, -R)
        except np.linalg.LinAlgError:
            print("Matrice singulière - arrêt")
            break
        
        y = y + delta
        
        if np.linalg.norm(delta, np.inf) < tol:
            print(f"Convergence en {iteration+1} itérations")
            break
    
    return x, y
```

Exemple d'utilisation pour l'équation de Bratu :

```python
# Définition de la non-linéarité
def bratu_f(x, y, yp):
    return -3.0 * np.exp(y)  # y'' = -3 exp(y)

conditions = [('D', 0, 0), ('D', -1, 0)]  # y(1)=0, y(-1)=0

x, y = solveur_nonlinear_generique(30, bratu_f, conditions)

plt.plot(x, y, 'b-', linewidth=2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Solution de Bratu avec solveur générique')
plt.grid(True)
plt.show()
```

---

## Section 4 : Analyse de convergence et validation

### 4.1. Convergence spectrale pour problèmes non-linéaires

La convergence spectrale se maintient pour les problèmes non-linéaires, comme le montre l'étude de l'équation de Bratu :

```python
def etude_convergence_bratu():
    """
    Étudie la convergence en fonction de N pour l'équation de Bratu
    """
    N_values = [5, 10, 15, 20, 25, 30, 35, 40]
    solutions = {}
    
    for N in N_values:
        x, y = bratu_nonlinear(N, 3.0)
        solutions[N] = (x, y)
    
    # Solution de référence (N maximal)
    x_ref, y_ref = bratu_nonlinear(50, 3.0)
    
    # Interpolation pour comparaison
    erreurs = []
    for N in N_values[:-1]:
        x_N, y_N = solutions[N]
        # Interpolation de y_N sur x_ref
        y_interp = np.interp(x_ref, x_N, y_N)
        erreur = np.max(np.abs(y_interp - y_ref))
        erreurs.append(erreur)
    
    # Graphique de convergence
    plt.figure(figsize=(10, 6))
    plt.semilogy(N_values[:-1], erreurs, 'bo-', linewidth=2, markersize=8)
    plt.xlabel('N (nombre de points)', fontsize=12)
    plt.ylabel('Erreur maximale', fontsize=12)
    plt.title('Convergence spectrale - Équation de Bratu non-linéaire', fontsize=14)
    plt.grid(True, alpha=0.3, which='both')
    plt.show()
    
    return N_values[:-1], erreurs
```

### 4.2. Comparaison linéaire vs non-linéaire

| Aspect | EDO linéaire | EDO non-linéaire |
|--------|--------------|------------------|
| **Système obtenu** | Linéaire $A\mathbf{y} = \mathbf{f}$ | Non-linéaire $F(\mathbf{y}) = 0$ |
| **Méthode de résolution** | Solveur direct unique | Itérations de Newton |
| **Sensibilité à l'initialisation** | Aucune | Critique |
| **Nombre de solutions** | Unique | Multiples (bifurcations) |
| **Coût de calcul** | $O(N^3)$ une fois | $O(N^3)$ par itération |
| **Difficulté d'implémentation** | Faible | Moyenne |

### 4.3. Avantages de la collocation pour les non-linéarités

La méthode de collocation présente des avantages décisifs pour les problèmes non-linéaires :

1. **Évaluation ponctuelle** : Les termes non-linéaires $f(y(x_j))$ sont simplement évalués aux nœuds
2. **Jacobienne accessible** : La matrice jacobienne a une structure simple
3. **Pas de transformation spectrale** : On reste dans l'espace physique
4. **Multiples solutions** : Différentes initialisations permettent de capturer différentes branches

---

## Conclusion

### Synthèse des résultats

La méthode de collocation de Tchebychev pour la résolution d'équations différentielles ordinaires **non-linéaires** du second ordre présente des avantages décisifs :

1. **Précision spectrale** maintenue même pour les non-linéarités fortes
2. **Simplicité d'implémentation** des termes non-linéaires
3. **Couplage naturel** avec la méthode de Newton
4. **Flexibilité** pour divers types de conditions aux limites

### Analyse comparative

| Méthode | Non-linéarités | Implémentation | Convergence | Domaine d'application |
|---------|----------------|----------------|-------------|----------------------|
| Galerkin | Très difficile | Complexe | Spectrale | Problèmes simples |
| Tau | Extrêmement difficile | Très complexe | Spectrale | Cas académiques |
| **Collocation** | **Très facile** | **Simple** | **Spectrale** | **Tous problèmes** |
| Diff. finies | Facile | Simple | Algébrique | Grands domaines |

### Perspectives

Les extensions naturelles de ce travail incluent :

1. **Problèmes avec bifurcations** : Continuation par paramètre pour suivre les branches de solutions
2. **Systèmes d'EDO non-linéaires** : Extension directe avec des inconnues vectorielles
3. **Problèmes raides** : Adaptation avec des méthodes de quasi-Newton
4. **Équations aux dérivées partielles non-linéaires** : Application aux EDP par produit tensoriel

La méthode de collocation de Tchebychev couplée à la méthode de Newton constitue ainsi un outil puissant et flexible pour la résolution numérique d'EDO non-linéaires du second ordre, alliant la précision spectrale à la simplicité d'implémentation des non-linéarités.