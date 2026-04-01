# Rapport de Recherche

## Résolution d'équations différentielles ordinaires du second ordre par méthodes de collocation de Tchebychev

---

## Introduction : Fondements des méthodes spectrales

Les méthodes spectrales constituent une classe avancée de techniques numériques pour la résolution d'équations différentielles. Contrairement aux approches locales telles que les différences finies ou les éléments finis, qui utilisent une information restreinte autour de chaque nœud, les méthodes spectrales reposent sur une représentation globale de la solution. Celle-ci est exprimée comme une combinaison linéaire de fonctions définies sur l’ensemble du domaine, ce qui permet d’obtenir une convergence extrêmement rapide pour les solutions régulières.

Le choix de la base de fonctions est déterminant. Les critères essentiels incluent la rapidité de convergence, la facilité de dérivation et d’intégration, l’orthogonalité, et la complétude. Pour les problèmes définis sur un intervalle borné et non périodiques, les polynômes de Tchebychev constituent un choix privilégié. Leur distribution des points de collocation, plus dense aux bords, minimise l’erreur d’interpolation et évite le phénomène de Runge, qui affecte les interpolations polynomiales sur des points équidistants.

---

## Section 1 : Méthodes de calcul des coefficients spectraux

La transformation d’une équation différentielle continue en un système algébrique discret repose sur l’approximation de la solution par une série tronquée :

$$u_N(x) = \sum_{k=0}^{N} \hat{u}_k \phi_k(x)$$

Le résidu associé est défini par $R(x; \hat{u}) = L u_N - f$, où $L$ est l’opérateur différentiel. Différentes stratégies permettent de minimiser ce résidu.

### 1.1. Méthode de Galerkin

Dans la méthode de Galerkin, on impose que le résidu soit orthogonal à chaque fonction de la base :

$$\langle R, \phi_j \rangle = \int_{\Omega} R(x) \phi_j(x) w(x) dx = 0, \quad j = 0, \dots, N$$

Cette approche exige que chaque fonction de base satisfasse individuellement les conditions aux limites. Elle conduit à des matrices symétriques et bien conditionnées, mais la construction d’une base adaptée peut s’avérer complexe, notamment pour des conditions aux limites non standards.

### 1.2. Méthode Tau de Lanczos

La méthode Tau est une extension de Galerkin pour les bases qui ne satisfont pas les conditions aux limites. Un terme correcteur est ajouté au résidu, et les dernières équations du système sont remplacées par les conditions aux limites. Cette approche travaille naturellement dans l’espace des coefficients, mais sa mise en œuvre est plus délicate et son extension aux problèmes non linéaires est difficile.

### 1.3. Méthode de collocation (pseudospectrale)

La méthode de collocation est la plus répandue en pratique. On exige que le résidu s’annule en un ensemble de points spécifiques, appelés nœuds de collocation :

$$R(x_j) = 0, \quad j = 0, \dots, N$$

Cela revient à utiliser des fonctions de test de type Dirac. Pour une équation du second ordre, la discrétisation conduit à un système algébrique impliquant les matrices de différenciation. Cette approche présente l’avantage majeur de traiter directement les non‑linéarités par évaluation ponctuelle, sans calcul d’intégrales complexes.

### 1.4. Synthèse comparative

| Méthode | Inconnues | Base de test | Structure matricielle | Non‑linéarités |
|---------|-----------|--------------|------------------------|----------------|
| Galerkin | Coefficients | $\phi_k$ (adaptée) | Bande | Difficile |
| Tau | Coefficients | $T_k$ | Presque bande | Très difficile |
| Collocation | Valeurs nodales | $\delta(x-x_j)$ | Pleine | Très facile |

---

## Section 2 : Propriétés des polynômes de Tchebychev

### 2.1. Définitions et relations fondamentales

Les polynômes de Tchebychev du premier genre $T_n(x)$ sont définis par la relation trigonométrique :

$$T_n(x) = \cos(n \arccos x), \quad x \in [-1, 1]$$

Ils satisfont la récurrence :

$$T_{n+1}(x) = 2x T_n(x) - T_{n-1}(x), \quad T_0(x) = 1, \quad T_1(x) = x$$

Le second genre $U_n(x)$ est défini par :

$$U_n(x) = \frac{\sin((n+1)\arccos x)}{\sin(\arccos x)}$$

### 2.2. Orthogonalité

Les polynômes de Tchebychev possèdent une propriété d’orthogonalité continue avec poids $1/\sqrt{1-x^2}$ :

$$\int_{-1}^{1} \frac{T_m(x) T_n(x)}{\sqrt{1-x^2}} dx = 
\begin{cases}
0 & m \neq n \\
\pi & m = n = 0 \\
\pi/2 & m = n \neq 0
\end{cases}$$

Ils vérifient également une orthogonalité discrète aux nœuds de collocation, propriété essentielle pour les calculs numériques.

### 2.3. Grilles de Tchebychev-Gauss-Lobatto

Les nœuds de collocation les plus couramment utilisés sont ceux de Gauss‑Lobatto :

$$x_j = \cos\left(\frac{\pi j}{N}\right), \quad j = 0, 1, \dots, N$$

Ces nœuds incluent les bornes du domaine et présentent une densité accrue aux extrémités, ce qui garantit la stabilité et la précision de l’interpolation.

---

## Section 3 : Opérateurs spectraux et résolution d’EDO

### 3.1. Matrice de différenciation de Tchebychev

La matrice de différenciation $D$ permet de calculer les dérivées aux nœuds de collocation. L’algorithme suivant construit cette matrice pour $N$ intervalles :

```python
import numpy as np

def cheb_matrix(N):
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0], c[-1] = 2, 2
    c = c * (-1)**np.arange(N + 1)
    X = np.tile(x, (N + 1, 1))
    dX = X - X.T
    D = np.outer(c, 1/c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x
```

La matrice de dérivée seconde s’obtient par $D^2 = D \cdot D$.

### 3.2. Résolution d’une équation linéaire

Considérons l’équation modèle :

$$y'' + y = 0, \quad y(-1) = 0, \quad y(1) = 1$$

La procédure de résolution est la suivante :

```python
N = 20
D, x = cheb_matrix(N)
D2 = D @ D
L = D2 + np.eye(N + 1)

L[0, :] = 0
L[0, 0] = 1
L[-1, :] = 0
L[-1, -1] = 1

F = np.zeros(N + 1)
F[0] = 1
F[-1] = 0

y = np.linalg.solve(L, F)
```

La solution peut être comparée à la solution analytique $y_{\text{exacte}} = \frac{\sin(x+1)}{\sin(2)}$.

### 3.3. Généralisation à une équation linéaire à coefficients variables

Pour une équation de la forme :

$$y''(x) + p(x) y'(x) + q(x) y(x) = f(x), \quad y(-1) = \alpha, \quad y(1) = \beta$$

on construit l’opérateur :

```python
p = p_fonction(x)
q = q_fonction(x)
f = f_fonction(x)

L = D2 + np.diag(p) @ D + np.diag(q)
L[0, :] = 0; L[0, 0] = 1
L[-1, :] = 0; L[-1, -1] = 1

F = f.copy()
F[0] = beta
F[-1] = alpha

y = np.linalg.solve(L, F)
```

### 3.4. Conditions de Neumann

Pour une condition de Neumann du type $y'(1) = \gamma$, on utilise la matrice $D$ :

```python
L[0, :] = D[0, :]
F[0] = gamma
```

### 3.5. Résolution d’équations non linéaires

La méthode de collocation se prête particulièrement bien aux problèmes non linéaires. Considérons l’équation de Bratu :

$$y''(x) + \lambda e^{y(x)} = 0, \quad y(-1) = y(1) = 0$$

Le système non linéaire $F(y) = 0$ est résolu par la méthode de Newton :

```python
def bratu_nonlinear(N, lambda_param):
    D, x = cheb_matrix(N)
    D2 = D @ D
    y = np.zeros(N + 1)
    tol = 1e-12
    for _ in range(50):
        F = D2 @ y + lambda_param * np.exp(y)
        J = D2.copy()
        for i in range(N + 1):
            J[i, i] += lambda_param * np.exp(y[i])
        J[0, :] = 0; J[0, 0] = 1; F[0] = y[0]
        J[-1, :] = 0; J[-1, -1] = 1; F[-1] = y[-1]
        delta = np.linalg.solve(J, -F)
        y += delta
        if np.linalg.norm(delta, np.inf) < tol:
            break
    return x, y
```

### 3.6. Équation de Duffing

L’équation de Duffing avec non linéarité cubique :

$$y''(x) + \omega^2 y(x) + \varepsilon y^3(x) = 0, \quad y(-1) = 0, \quad y(1) = 1$$

est résolue de manière analogue, avec une jacobienne tenant compte de la dérivée $3\varepsilon y^2$.

### 3.7. Problème mixte avec conditions de Neumann

Pour une équation du type :

$$y''(x) + y'(x) + y^2(x) = \cos(\pi x), \quad y'(-1) = 0, \quad y(1) = 0$$

la jacobienne inclut la contribution de la dérivée première et de la non linéarité quadratique. La condition de Neumann est imposée en remplaçant la ligne correspondante par celle de la matrice $D$.

---

## Section 4 : Analyse de convergence

### 4.1. Convergence spectrale

La méthode de collocation de Tchebychev converge exponentiellement pour les solutions régulières. L’erreur maximale décroît plus vite que toute puissance polynomiale de $1/N$. Pour l’équation de Bratu, une étude systématique montre qu’une précision machine est atteinte pour $N$ de l’ordre de 30 à 40.

### 4.2. Comparaison avec les différences finies

| Méthode | Points pour $\varepsilon = 10^{-10}$ | Complexité |
|---------|--------------------------------------|------------|
| Différences finies (ordre 2) | $\sim 10^6$ | $O(N)$ |
| Différences finies (ordre 4) | $\sim 10^4$ | $O(N)$ |
| Collocation Tchebychev | $\sim 30$ | $O(N^3)$ (direct) |

Bien que la collocation utilise des matrices pleines, le faible nombre de points nécessaire la rend souvent plus efficace globalement pour les problèmes nécessitant une haute précision.

---

## Conclusion

La méthode de collocation de Tchebychev constitue une approche puissante pour la résolution d’équations différentielles ordinaires du second ordre. Ses principaux atouts sont les suivants :

- **Convergence spectrale** : une précision exceptionnelle avec un nombre réduit de points.
- **Simplicité d’implémentation** : les non‑linéarités sont traitées par évaluation ponctuelle.
- **Flexibilité** : elle s’adapte facilement aux coefficients variables et aux conditions aux limites mixtes.
- **Stabilité** : l’utilisation des nœuds de Tchebychev élimine le phénomène de Runge.

La comparaison avec les méthodes de Galerkin et Tau montre que la collocation est la plus adaptée aux problèmes complexes, notamment non linéaires, au prix de matrices pleines et d’un conditionnement qui se dégrade pour de très grands $N$.

Parmi les perspectives, on peut citer l’extension aux systèmes d’équations, aux équations aux dérivées partielles par produit tensoriel, l’utilisation de transformées rapides pour réduire la complexité, et la continuation pour suivre des branches de solutions en présence de bifurcations.

---

## Références

1. Trefethen, L. N. (2000). *Spectral Methods in MATLAB*. SIAM.
2. Boyd, J. P. (2001). *Chebyshev and Fourier Spectral Methods*. Dover Publications.
3. Canuto, C., Hussaini, M. Y., Quarteroni, A., & Zang, T. A. (2006). *Spectral Methods: Fundamentals in Single Domains*. Springer.
4. Gottlieb, D., & Orszag, S. A. (1977). *Numerical Analysis of Spectral Methods: Theory and Applications*. SIAM.