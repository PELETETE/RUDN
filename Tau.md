## Méthode Tau de Lanczos

### Principe général

La méthode Tau, introduite par Lanczos, est une extension de la méthode de Galerkin qui permet d'utiliser une base de fonctions ne satisfaisant pas individuellement les conditions aux limites. L'idée fondamentale est d'ajouter un terme correcteur polynomial (le "tau") au résidu pour compenser le non-respect des conditions aux limites.

Soit l'équation différentielle linéaire d'ordre 2 :

$$L u(x) = f(x), \quad x \in [-1, 1]$$

avec des conditions aux limites :

$$u(-1) = \alpha, \quad u(1) = \beta$$

### Développement en série de Chebyshev

On cherche une solution approchée sous forme d'une série tronquée de polynômes de Chebyshev :

$$u_N(x) = \sum_{k=0}^{N} \hat{u}_k T_k(x)$$

où $T_k(x)$ sont les polynômes de Chebyshev de première espèce, définis par $T_k(x) = \cos(k \arccos x)$.

### Formulation du résidu

Le résidu de l'équation différentielle est :

$$R(x) = L u_N(x) - f(x)$$

L'idée de Lanczos est d'ajouter un terme correcteur polynomial de degré $N$ ou $N+1$ pour que les conditions aux limites puissent être satisfaites. On écrit :

$$L u_N(x) - f(x) = \tau_1 \psi_1(x) + \tau_2 \psi_2(x) + \cdots$$

où $\psi_i(x)$ sont des polynômes de degré élevé (typiquement $T_N(x)$ et $T_{N-1}(x)$).

### Projection de Galerkin

On projette l'équation modifiée sur les $N+1$ polynômes de base $T_j(x)$ avec $j = 0, \dots, N$ :

$$\int_{-1}^1 [L u_N(x) - f(x)] T_j(x) \frac{dx}{\sqrt{1-x^2}} = \int_{-1}^1 \left[\sum_{i=1}^{m} \tau_i \psi_i(x)\right] T_j(x) \frac{dx}{\sqrt{1-x^2}}$$

En utilisant l'orthogonalité des polynômes de Chebyshev :

$$\int_{-1}^1 T_j(x) T_k(x) \frac{dx}{\sqrt{1-x^2}} = 
\begin{cases}
0 & j \neq k \\
\pi & j = k = 0 \\
\pi/2 & j = k \neq 0
\end{cases}$$

### Système d'équations

Pour un opérateur linéaire $L$ du second ordre, l'application de $L$ aux polynômes de Chebyshev s'exprime par des relations de récurrence. On obtient un système linéaire de $N+1$ équations pour les coefficients $\hat{u}_k$ et les paramètres $\tau_i$ :

$$\sum_{k=0}^{N} \hat{u}_k \langle L T_k, T_j \rangle - \langle f, T_j \rangle = \sum_{i=1}^{m} \tau_i \langle \psi_i, T_j \rangle, \quad j = 0, \dots, N$$

où $\langle \cdot, \cdot \rangle$ désigne le produit scalaire avec poids $w(x) = 1/\sqrt{1-x^2}$.

### Ajout des conditions aux limites

Les $N+1$ équations de projection sont modifiées : on remplace les $m$ dernières équations (celles de plus haut degré) par les conditions aux limites.

Pour un problème du second ordre avec deux conditions aux limites, on utilise généralement $m=2$ termes de correction. Les conditions aux limites s'écrivent :

$$u_N(-1) = \sum_{k=0}^{N} \hat{u}_k T_k(-1) = \alpha$$

$$u_N(1) = \sum_{k=0}^{N} \hat{u}_k T_k(1) = \beta$$

En utilisant $T_k(1) = 1$ et $T_k(-1) = (-1)^k$, on obtient :

$$\sum_{k=0}^{N} \hat{u}_k = \beta$$

$$\sum_{k=0}^{N} (-1)^k \hat{u}_k = \alpha$$

### Structure du système final

Le système complet s'écrit sous forme matricielle :

$$\begin{bmatrix}
A_{0,0} & A_{0,1} & \cdots & A_{0,N} & -\langle \psi_1, T_0 \rangle & -\langle \psi_2, T_0 \rangle \\
A_{1,0} & A_{1,1} & \cdots & A_{1,N} & -\langle \psi_1, T_1 \rangle & -\langle \psi_2, T_1 \rangle \\
\vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
A_{N-2,0} & A_{N-2,1} & \cdots & A_{N-2,N} & -\langle \psi_1, T_{N-2} \rangle & -\langle \psi_2, T_{N-2} \rangle \\
1 & 1 & \cdots & 1 & 0 & 0 \\
(-1)^0 & (-1)^1 & \cdots & (-1)^N & 0 & 0
\end{bmatrix}
\begin{bmatrix}
\hat{u}_0 \\ \hat{u}_1 \\ \vdots \\ \hat{u}_N \\ \tau_1 \\ \tau_2
\end{bmatrix}
=
\begin{bmatrix}
\langle f, T_0 \rangle \\
\langle f, T_1 \rangle \\
\vdots \\
\langle f, T_{N-2} \rangle \\
\beta \\
\alpha
\end{bmatrix}$$

où $A_{j,k} = \langle L T_k, T_j \rangle$.

### Choix des polynômes de correction

Les polynômes $\psi_i$ sont typiquement choisis comme les polynômes de Chebyshev de plus haut degré :

$$\psi_1(x) = T_N(x), \quad \psi_2(x) = T_{N-1}(x)$$

Avec ce choix, les produits scalaires $\langle \psi_i, T_j \rangle$ sont nuls pour $j < N-1$ en raison de l'orthogonalité, ce qui simplifie la structure du système.

### Équivalence avec la méthode de collocation

Pour certains opérateurs, la méthode Tau peut être reformulée comme une méthode de collocation aux nœuds de Gauss-Lobatto, offrant une interprétation alternative dans l'espace physique.

### Avantages et inconvénients

| Aspect | Caractéristique |
|--------|-----------------|
| **Avantages** | Travail direct dans l'espace spectral, pas de construction de base adaptée, matrice presque bande |
| **Inconvénients** | Gestion complexe des non-linéarités, nécessite le calcul d'intégrales pour les produits scalaires |