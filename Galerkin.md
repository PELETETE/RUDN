## Méthode de Galerkin pour la résolution d'EDO du second ordre par polynômes de Tchebychev

---

### 1. Principe général

La méthode de Galerkin est une méthode de projection où l'on cherche une solution approchée $u_N(x)$ sous forme d'une combinaison linéaire de fonctions de base $\phi_k(x)$, et l'on force le résidu à être orthogonal à l'espace de ces fonctions de base.

Soit l'équation différentielle générale du second ordre :

$$
\mathcal{L}u = f(x), \quad x \in [-1, 1]
$$

avec conditions aux limites :

$$
u(-1) = \alpha, \quad u(1) = \beta
$$

On approxime la solution par :

$$
u_N(x) = \sum_{k=0}^{N} \hat{u}_k \phi_k(x)
$$

où $\phi_k(x)$ sont des fonctions de base choisies.

Le résidu est défini par :

$$
R_N(x) = \mathcal{L}u_N(x) - f(x)
$$

La méthode de Galerkin impose l'orthogonalité du résidu à chaque fonction de base :

$$
\int_{-1}^{1} R_N(x) \phi_j(x) w(x) \, dx = 0, \quad j = 0, 1, \dots, N
$$

où $w(x)$ est une fonction de poids (pour les polynômes de Tchebychev, $w(x) = 1/\sqrt{1-x^2}$).

---

### 2. Base adaptée aux conditions aux limites

Pour les polynômes de Tchebychev, un problème crucial se pose : les polynômes $T_k(x)$ ne satisfont pas les conditions aux limites $T_k(\pm 1) \neq 0$ pour $k$ pair.

Il est donc nécessaire de construire une **base adaptée** qui satisfait automatiquement les conditions aux limites homogènes $u(\pm 1) = 0$.

#### 2.1. Base pour conditions de Dirichlet homogènes

Pour $u(-1) = u(1) = 0$, on utilise :

$$
\phi_k(x) = T_{k+2}(x) - T_k(x), \quad k = 0, 1, \dots, N-2
$$

**Propriété** : 
$$
\phi_k(1) = T_{k+2}(1) - T_k(1) = 1 - 1 = 0
$$
$$
\phi_k(-1) = T_{k+2}(-1) - T_k(-1) = (-1)^{k+2} - (-1)^k = (-1)^k - (-1)^k = 0
$$

La solution approchée s'écrit alors :

$$
u_N(x) = \sum_{k=0}^{N-2} \hat{u}_k \left[ T_{k+2}(x) - T_k(x) \right]
$$

#### 2.2. Base pour conditions de Dirichlet non-homogènes

Pour $u(-1) = \alpha$, $u(1) = \beta$, on décompose :

$$
u_N(x) = u_0(x) + \tilde{u}_N(x)
$$

où $u_0(x)$ est une fonction satisfaisant les conditions aux limites, par exemple la droite :

$$
u_0(x) = \frac{\alpha + \beta}{2} + \frac{\beta - \alpha}{2}x
$$

et $\tilde{u}_N(x)$ est la solution de l'équation homogène avec $\tilde{u}_N(\pm 1) = 0$, développée sur la base adaptée.

---

### 3. Formulation variationnelle

#### 3.1. Opérateur différentiel général

Considérons l'EDO :

$$
u''(x) + p(x)u'(x) + q(x)u(x) = f(x)
$$

En remplaçant $u$ par $u_N = \sum_{k=0}^{M} \hat{u}_k \phi_k(x)$ (avec $M = N-2$), on obtient :

$$
\sum_{k=0}^{M} \hat{u}_k \left[ \phi_k''(x) + p(x)\phi_k'(x) + q(x)\phi_k(x) \right] = f(x) + R_N(x)
$$

#### 3.2. Projection de Galerkin

Les équations de Galerkin s'écrivent :

$$
\sum_{k=0}^{M} \hat{u}_k \int_{-1}^{1} \left[ \phi_k''(x) + p(x)\phi_k'(x) + q(x)\phi_k(x) \right] \phi_j(x) w(x) \, dx = \int_{-1}^{1} f(x) \phi_j(x) w(x) \, dx
$$

pour $j = 0, 1, \dots, M$.

Sous forme matricielle :

$$
\mathbf{A} \hat{\mathbf{u}} = \mathbf{b}
$$

avec :

$$
A_{jk} = \int_{-1}^{1} \left[ \phi_k''(x) + p(x)\phi_k'(x) + q(x)\phi_k(x) \right] \phi_j(x) w(x) \, dx
$$

$$
b_j = \int_{-1}^{1} f(x) \phi_j(x) w(x) \, dx
$$

---

### 4. Intégration par parties pour la symétrie

Pour obtenir des matrices symétriques (propriété souhaitable), on intègre par parties le terme du second ordre :

$$
\int_{-1}^{1} \phi_k''(x) \phi_j(x) w(x) \, dx = \left[ \phi_k'(x) \phi_j(x) w(x) \right]_{-1}^{1} - \int_{-1}^{1} \phi_k'(x) \left[ \phi_j'(x) w(x) + \phi_j(x) w'(x) \right] dx
$$

Pour le poids de Tchebychev $w(x) = 1/\sqrt{1-x^2}$, les termes de bord s'annulent car $w(\pm 1) \to \infty$ mais avec un comportement tel que le produit reste fini.

La forme faible devient :

$$
\int_{-1}^{1} \phi_k'(x) \phi_j'(x) w(x) \, dx - \int_{-1}^{1} \phi_k'(x) \phi_j(x) w'(x) \, dx + \int_{-1}^{1} \left[ p(x)\phi_k'(x) + q(x)\phi_k(x) \right] \phi_j(x) w(x) \, dx
$$

---

### 5. Calcul pratique des intégrales

#### 5.1. Intégration par quadrature de Gauss-Chebyshev

Les intégrales sont calculées par quadrature de Gauss-Chebyshev. Pour le poids $w(x) = 1/\sqrt{1-x^2}$, les points et poids de quadrature sont :

**Points de Gauss-Chebyshev** (racines de $T_{M+1}(x)$) :

$$
x_i = \cos\left( \frac{(2i+1)\pi}{2(M+1)} \right), \quad i = 0, \dots, M
$$

**Poids** :

$$
w_i = \frac{\pi}{M+1}
$$

L'intégrale s'approxime par :

$$
\int_{-1}^{1} g(x) \frac{dx}{\sqrt{1-x^2}} \approx \sum_{i=0}^{M} g(x_i) \frac{\pi}{M+1}
$$

#### 5.2. Cas particulier : opérateur $u'' + \lambda u$

Pour l'équation modèle $u'' + \lambda u = f$ avec $u(\pm 1) = 0$, la matrice de Galerkin devient :

$$
A_{jk} = \int_{-1}^{1} \phi_k''(x) \phi_j(x) w(x) \, dx + \lambda \int_{-1}^{1} \phi_k(x) \phi_j(x) w(x) \, dx
$$

Après intégration par parties :

$$
A_{jk} = -\int_{-1}^{1} \phi_k'(x) \phi_j'(x) w(x) \, dx + \lambda \int_{-1}^{1} \phi_k(x) \phi_j(x) w(x) \, dx
$$

La matrice obtenue est symétrique.

---

### 6. Relation entre les coefficients et les valeurs nodales

Il est parfois utile de passer des coefficients spectraux $\hat{u}_k$ aux valeurs nodales $u(x_i)$. La matrice de passage est donnée par :

$$
u(x_i) = \sum_{k=0}^{M} \hat{u}_k \phi_k(x_i)
$$

ou sous forme matricielle :

$$
\mathbf{u} = \mathbf{\Phi} \hat{\mathbf{u}}
$$

où $\Phi_{ik} = \phi_k(x_i)$.

La transformation inverse (si la matrice est inversible) permet de reconstruire la solution en tout point.

---

### 7. Algorithme de résolution

1. **Choix de la base** : Construire $\phi_k(x) = T_{k+2}(x) - T_k(x)$ pour $k = 0, \dots, N-2$
2. **Discrétisation** : Choisir $M+1$ points de quadrature de Gauss-Chebyshev
3. **Construction de la matrice** : Calculer $A_{jk}$ par quadrature
   $$
   A_{jk} = \sum_{i=0}^{M} \left[ \phi_k''(x_i) + p(x_i)\phi_k'(x_i) + q(x_i)\phi_k(x_i) \right] \phi_j(x_i) w_i
   $$
4. **Construction du second membre** : Calculer $b_j = \sum_{i=0}^{M} f(x_i) \phi_j(x_i) w_i$
5. **Résolution** : Résoudre $\mathbf{A} \hat{\mathbf{u}} = \mathbf{b}$ pour obtenir les coefficients $\hat{u}_k$
6. **Reconstruction** : Calculer $u_N(x)$ en tout point via $u_N(x) = \sum_{k=0}^{M} \hat{u}_k \phi_k(x)$

---

### 8. Avantages et inconvénients

#### Avantages
- Matrice symétrique pour les opérateurs auto-adjoints
- Meilleur conditionnement que la collocation
- Convergence spectrale garantie
- Base adaptée aux conditions aux limites

#### Inconvénients
- Nécessité de construire une base adaptée (différente selon les conditions aux limites)
- Calculs d'intégrales potentiellement coûteux
- Gestion complexe des non-linéarités (nécessite des produits de convolution)
- Difficulté à traiter les coefficients variables rapidement

---

### 9. Exemple : $u'' + u = 0$, $u(\pm 1) = 0$

Pour $p(x) = 0$, $q(x) = 1$, $f(x) = 0$, la matrice de Galerkin est :

$$
A_{jk} = -\int_{-1}^{1} \phi_k'(x) \phi_j'(x) \frac{dx}{\sqrt{1-x^2}} + \int_{-1}^{1} \phi_k(x) \phi_j(x) \frac{dx}{\sqrt{1-x^2}}
$$

Les intégrales se calculent analytiquement en utilisant les propriétés d'orthogonalité des polynômes de Tchebychev. On obtient une matrice tridiagonale dont les coefficients sont donnés par des formules fermées.