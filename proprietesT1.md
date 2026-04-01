```markdown
# Polynômes de Tchebychev de première espèce : propriétés et implémentation

---

## 1. Définition fondamentale

Les polynômes de Tchebychev de première espèce $T_n(x)$ sont définis sur l'intervalle $[-1, 1]$ par la relation trigonométrique :

$$T_n(x) = \cos(n \arccos x), \quad n \in \mathbb{N}, \quad x \in [-1, 1]$$

Cette définition est valable pour tout $n \geq 0$. Pour $n = 0$ et $n = 1$, on a :

$$T_0(x) = \cos(0) = 1$$
$$T_1(x) = \cos(\arccos x) = x$$

---

## 2. Relation de récurrence

Les polynômes de Tchebychev satisfont la relation de récurrence à trois termes :

$$T_{n+1}(x) = 2x T_n(x) - T_{n-1}(x), \quad n \geq 1$$

Cette relation permet de générer tous les polynômes de manière récursive.

---

## 3. Forme explicite

Pour $n \geq 0$, on peut écrire :

$$T_n(x) = \frac{n}{2} \sum_{k=0}^{\lfloor n/2 \rfloor} (-1)^k \frac{(n-k-1)!}{k!(n-2k)!} (2x)^{n-2k}$$

ou encore :

$$T_n(x) = \frac{1}{2} \left[ (x + \sqrt{x^2 - 1})^n + (x - \sqrt{x^2 - 1})^n \right]$$

---

## 4. Polynômes pour les premiers degrés

Les premiers polynômes de Tchebychev sont :

$$T_0(x) = 1$$
$$T_1(x) = x$$
$$T_2(x) = 2x^2 - 1$$
$$T_3(x) = 4x^3 - 3x$$
$$T_4(x) = 8x^4 - 8x^2 + 1$$
$$T_5(x) = 16x^5 - 20x^3 + 5x$$
$$T_6(x) = 32x^6 - 48x^4 + 18x^2 - 1$$

---

## 5. Orthogonalité

Les polynômes de Tchebychev sont orthogonaux sur $[-1, 1]$ par rapport à la fonction poids $w(x) = 1/\sqrt{1-x^2}$ :

$$\int_{-1}^{1} T_m(x) T_n(x) \frac{dx}{\sqrt{1-x^2}} = 
\begin{cases}
0 & m \neq n \\
\pi & m = n = 0 \\
\frac{\pi}{2} & m = n \neq 0
\end{cases}$$

---

## 6. Orthogonalité discrète

Aux nœuds de Tchebychev-Gauss-Lobatto $x_j = \cos(j\pi/N)$, $j = 0, \dots, N$, on a la propriété d'orthogonalité discrète :

$$\sum_{j=0}^{N} T_m(x_j) T_n(x_j) = 
\begin{cases}
0 & m \neq n \\
N+1 & m = n = 0 \\
\frac{N+1}{2} & 0 < m = n < N \\
N & m = n = N
\end{cases}$$

où les sommes sont pondérées par des facteurs $c_j = 2$ pour $j=0,N$ et $c_j = 1$ sinon.

---

## 7. Propriété minimax

Parmi tous les polynômes de degré $n$ avec coefficient dominant $2^{n-1}$, $T_n(x)$ minimise la norme infinie sur $[-1, 1]$ :

$$\max_{x \in [-1,1]} |T_n(x)| = 1$$

Les extrema de $T_n(x)$ sont atteints aux points $x_k = \cos(k\pi/n)$ avec $T_n(x_k) = (-1)^k$.

---

## 8. Dérivées

La dérivée première s'exprime par :

$$T_n'(x) = n U_{n-1}(x)$$

où $U_{n-1}$ est le polynôme de Tchebychev de seconde espèce.

Plus généralement :

$$\frac{d^k}{dx^k} T_n(x) = 2^{k-1} n \prod_{j=0}^{k-1} (n^2 - j^2) \frac{\sin((n-k)\arccos x)}{\sin^k(\arccos x)} + \dots$$

Pour les nœuds de collocation, on utilise les matrices de différenciation.

---

## 9. Racines et nœuds

Les racines de $T_n(x)$ sont données par :

$$x_k = \cos\left(\frac{2k-1}{2n}\pi\right), \quad k = 1, \dots, n$$

Les extrema (nœuds de Gauss-Lobatto) sont :

$$x_k = \cos\left(\frac{k\pi}{n}\right), \quad k = 0, \dots, n$$

---

## 10. Expansion en série de Tchebychev

Toute fonction $f(x)$ définie sur $[-1, 1]$ peut être approchée par une série tronquée :

$$f(x) \approx \sum_{k=0}^{N} a_k T_k(x)$$

Les coefficients sont donnés par :

$$a_0 = \frac{1}{\pi} \int_{-1}^{1} f(x) \frac{dx}{\sqrt{1-x^2}}$$
$$a_k = \frac{2}{\pi} \int_{-1}^{1} f(x) T_k(x) \frac{dx}{\sqrt{1-x^2}}, \quad k \geq 1$$

---

## 11. Approximation polynomiale

L'erreur d'interpolation polynomiale aux nœuds de Tchebychev est minimisée :

$$\|f - p_N\|_\infty \leq \frac{2}{\pi} \log(N+1) \inf_{q \in \mathcal{P}_N} \|f - q\|_\infty$$

où $\mathcal{P}_N$ est l'espace des polynômes de degré $\leq N$.

---

## 12. Transformation de domaine

Pour un problème sur $[a, b]$, on utilise la transformation affine :

$$x = \frac{2}{b-a} t - \frac{a+b}{b-a}, \quad t \in [a, b]$$

soit :

$$t = \frac{b-a}{2} x + \frac{a+b}{2}$$

Les polynômes de Tchebychev s'adaptent alors à l'intervalle $[a, b]$ via :

$$\tilde{T}_n(t) = T_n\left(\frac{2t - (a+b)}{b-a}\right)$$

---

## 13. Code : Génération et propriétés

```python
import numpy as np
import matplotlib.pyplot as plt

def chebyshev_poly(n, x):
    """Calcule T_n(x) pour n donné et un vecteur x"""
    if n == 0:
        return np.ones_like(x)
    if n == 1:
        return x
    T_n_2 = np.ones_like(x)
    T_n_1 = x.copy()
    for k in range(2, n + 1):
        T_n = 2 * x * T_n_1 - T_n_2
        T_n_2, T_n_1 = T_n_1, T_n
    return T_n_1

def chebyshev_roots(n):
    """Racines de T_n(x)"""
    k = np.arange(1, n + 1)
    return np.cos((2 * k - 1) * np.pi / (2 * n))

def chebyshev_extrema(n):
    """Extrema de T_n(x) (nœuds de Gauss-Lobatto)"""
    return np.cos(np.arange(n + 1) * np.pi / n)

def chebyshev_coeffs(f, N):
    """Coefficients de Tchebychev par quadrature"""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    fx = f(x)
    c = np.ones(N + 1)
    c[0], c[-1] = 2, 2
    a = np.zeros(N + 1)
    for k in range(N + 1):
        Tk = np.cos(k * np.arccos(x))
        a[k] = (2 / (N * c[k])) * np.sum(fx * Tk / c)
    a[0] /= 2
    return a

def chebyshev_eval(a, x):
    """Évaluation de la série de Tchebychev"""
    N = len(a) - 1
    if N == 0:
        return a[0] * np.ones_like(x)
    T_n_2 = np.ones_like(x)
    T_n_1 = x.copy()
    result = a[0] * T_n_2 + a[1] * T_n_1
    for n in range(2, N + 1):
        T_n = 2 * x * T_n_1 - T_n_2
        result += a[n] * T_n
        T_n_2, T_n_1 = T_n_1, T_n
    return result

def chebyshev_matrix(N):
    """Matrice de différenciation de Tchebychev"""
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

def plot_chebyshev_polynomials(n_max, N_points=500):
    """Trace les polynômes de Tchebychev T_0 à T_n_max"""
    x = np.linspace(-1, 1, N_points)
    plt.figure(figsize=(10, 6))
    for n in range(n_max + 1):
        y = chebyshev_poly(n, x)
        plt.plot(x, y, label=f'T_{n}(x)')
    plt.xlabel('x')
    plt.ylabel('T_n(x)')
    plt.title('Polynômes de Tchebychev de première espèce')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
```

---

## 14. Code : Orthogonalité et quadrature

```python
def check_orthogonality(N, n, m):
    """Vérifie l'orthogonalité continue de T_n et T_m"""
    x = np.linspace(-1, 1, 10000)
    dx = x[1] - x[0]
    weight = 1 / np.sqrt(1 - x**2)
    Tn = chebyshev_poly(n, x)
    Tm = chebyshev_poly(m, x)
    integral = np.sum(Tn * Tm * weight) * dx
    return integral

def chebyshev_quadrature(f, N):
    """Quadrature de Gauss-Tchebychev"""
    x = np.cos(np.pi * (2 * np.arange(N) + 1) / (2 * N))
    return (np.pi / N) * np.sum(f(x))

def discrete_orthogonality(N, n, m):
    """Vérifie l'orthogonalité discrète aux nœuds de Gauss-Lobatto"""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0], c[-1] = 2, 2
    Tn = np.cos(n * np.arccos(x))
    Tm = np.cos(m * np.arccos(x))
    return np.sum(Tn * Tm / c)
```

---

## 15. Code : Approximation de fonctions

```python
def approximate_function(f, N, x_eval=None):
    """Approxime f par une série de Tchebychev"""
    if x_eval is None:
        x_eval = np.linspace(-1, 1, 500)
    a = chebyshev_coeffs(f, N)
    return chebyshev_eval(a, x_eval)

# Exemple : approximation de exp(x)
def test_approximation():
    f = np.exp
    N_values = [5, 10, 15, 20]
    x = np.linspace(-1, 1, 500)
    y_exact = f(x)
    
    plt.figure(figsize=(12, 8))
    for N in N_values:
        y_approx = approximate_function(f, N, x)
        plt.subplot(2, 2, N_values.index(N) + 1)
        plt.plot(x, y_exact, 'b-', label='Exact')
        plt.plot(x, y_approx, 'r--', label=f'N={N}')
        plt.xlabel('x')
        plt.ylabel('exp(x)')
        plt.title(f'Approximation avec N={N}')
        plt.legend()
        plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
```

---

## 16. Code : Matrice de différenciation

```python
def compute_derivative_matrix(N):
    """Calcule la matrice de différenciation et vérifie les dérivées"""
    D, x = chebyshev_matrix(N)
    D2 = D @ D
    
    # Test sur une fonction test
    f = np.exp
    y = f(x)
    y_prime_exact = np.exp(x)
    y_prime_num = D @ y
    y_second_exact = np.exp(x)
    y_second_num = D2 @ y
    
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(x, y_prime_exact, 'b-', label='Exact')
    plt.plot(x, y_prime_num, 'ro', label='Numérique')
    plt.xlabel('x')
    plt.ylabel("y'(x)")
    plt.title('Dérivée première')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.plot(x, y_second_exact, 'b-', label='Exact')
    plt.plot(x, y_second_num, 'ro', label='Numérique')
    plt.xlabel('x')
    plt.ylabel("y''(x)")
    plt.title('Dérivée seconde')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    
    erreur_1 = np.max(np.abs(y_prime_num - y_prime_exact))
    erreur_2 = np.max(np.abs(y_second_num - y_second_exact))
    print(f"Erreur max dérivée première : {erreur_1:.2e}")
    print(f"Erreur max dérivée seconde : {erreur_2:.2e}")
    
    return D, x
```

---

## 17. Code : Convergence spectrale

```python
def study_convergence(f, N_max, x_test=None):
    """Étude de la convergence de l'approximation"""
    if x_test is None:
        x_test = np.linspace(-1, 1, 100)
    y_exact = f(x_test)
    N_values = range(2, N_max + 1, 2)
    erreurs = []
    
    for N in N_values:
        a = chebyshev_coeffs(f, N)
        y_approx = chebyshev_eval(a, x_test)
        erreur = np.max(np.abs(y_approx - y_exact))
        erreurs.append(erreur)
    
    plt.figure(figsize=(10, 6))
    plt.semilogy(N_values, erreurs, 'bo-')
    plt.xlabel('N (degré)')
    plt.ylabel('Erreur maximale')
    plt.title('Convergence spectrale de l\'approximation de Tchebychev')
    plt.grid(True, alpha=0.3)
    plt.show()
    
    return N_values, erreurs

# Test sur une fonction analytique
study_convergence(np.exp, 40)
```

---

## 18. Code : Interpolation aux nœuds de Tchebychev

```python
def compare_interpolation(f, N):
    """Compare l'interpolation aux nœuds de Tchebychev vs équidistants"""
    x_cheb = np.cos(np.pi * np.arange(N + 1) / N)
    x_eq = np.linspace(-1, 1, N + 1)
    x_fine = np.linspace(-1, 1, 500)
    
    # Interpolation polynomiale
    from numpy.polynomial import Polynomial
    p_cheb = Polynomial.fit(x_cheb, f(x_cheb), N)
    p_eq = Polynomial.fit(x_eq, f(x_eq), N)
    
    y_cheb = p_cheb(x_fine)
    y_eq = p_eq(x_fine)
    y_exact = f(x_fine)
    
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(x_fine, y_exact, 'b-', label='Exact')
    plt.plot(x_fine, y_cheb, 'r--', label='Tchebychev')
    plt.plot(x_cheb, f(x_cheb), 'ro', markersize=4)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title(f'Interpolation aux nœuds de Tchebychev (N={N})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.plot(x_fine, y_exact, 'b-', label='Exact')
    plt.plot(x_fine, y_eq, 'g--', label='Équidistants')
    plt.plot(x_eq, f(x_eq), 'go', markersize=4)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title(f'Interpolation aux nœuds équidistants (N={N})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
```

---

## 19. Code : Transformation d'intervalle

```python
def map_to_interval(a, b):
    """Crée des fonctions de Tchebychev adaptées à [a, b]"""
    def transform(x):
        return (2 * x - (a + b)) / (b - a)
    
    def inverse_transform(t):
        return (b - a) * t / 2 + (a + b) / 2
    
    def chebyshev_on_interval(n, x):
        t = transform(x)
        return chebyshev_poly(n, t)
    
    return chebyshev_on_interval, transform, inverse_transform

# Exemple : T_3 adapté à [0, 1]
def test_interval_mapping():
    a, b = 0, 1
    cheb_on_interval, transform, _ = map_to_interval(a, b)
    x = np.linspace(a, b, 200)
    y = cheb_on_interval(3, x)
    
    plt.figure(figsize=(8, 5))
    plt.plot(x, y, 'b-', linewidth=2)
    plt.xlabel('x')
    plt.ylabel('T_3(x) adapté')
    plt.title(f'Polynôme de Tchebychev T_3 sur [{a}, {b}]')
    plt.grid(True, alpha=0.3)
    plt.show()
```

---

## 20. Code : Applications supplémentaires

```python
def generate_chebyshev_table(n_max):
    """Génère un tableau des polynômes T_n(x)"""
    from sympy import symbols, chebyshevt, expand
    x = symbols('x')
    print("Polynômes de Tchebychev de première espèce :")
    print("-" * 40)
    for n in range(n_max + 1):
        Tn = expand(chebyshevt(n, x))
        print(f"T_{n}(x) = {Tn}")

def plot_roots_extrema(n):
    """Trace les racines et extrema de T_n(x)"""
    x_fine = np.linspace(-1, 1, 500)
    y = chebyshev_poly(n, x_fine)
    roots = chebyshev_roots(n)
    extrema = chebyshev_extrema(n)
    y_extrema = chebyshev_poly(n, extrema)
    
    plt.figure(figsize=(10, 6))
    plt.plot(x_fine, y, 'b-', linewidth=1.5, label=f'T_{n}(x)')
    plt.plot(roots, np.zeros_like(roots), 'ro', label='Racines')
    plt.plot(extrema, y_extrema, 'gs', label='Extrema')
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    plt.xlabel('x')
    plt.ylabel(f'T_{n}(x)')
    plt.title(f'Racines et extrema de T_{n}(x)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
```

---

## 21. Code : Exemple complet d'application

```python
def full_demo():
    """Démonstration complète des propriétés des polynômes de Tchebychev"""
    print("=== Polynômes de Tchebychev de première espèce ===\n")
    
    # 1. Génération des polynômes
    print("1. Génération des premiers polynômes :")
    x = np.linspace(-1, 1, 100)
    plt.figure(figsize=(10, 5))
    for n in range(5):
        plt.plot(x, chebyshev_poly(n, x), label=f'T_{n}(x)')
    plt.xlabel('x')
    plt.ylabel('T_n(x)')
    plt.title('Polynômes de Tchebychev')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
    
    # 2. Racines et extrema
    print("\n2. Racines et extrema de T_8(x) :")
    plot_roots_extrema(8)
    
    # 3. Approximation d'une fonction
    print("\n3. Approximation de exp(x) par série de Tchebychev :")
    test_approximation()
    
    # 4. Matrice de différenciation
    print("\n4. Matrice de différenciation et calcul des dérivées :")
    compute_derivative_matrix(15)
    
    # 5. Convergence
    print("\n5. Étude de convergence :")
    study_convergence(np.exp, 30)
    
    # 6. Interpolation comparée
    print("\n6. Comparaison Tchebychev vs équidistants :")
    def runge(x):
        return 1 / (1 + 25 * x**2)
    compare_interpolation(runge, 15)

if __name__ == "__main__":
    full_demo()
```

---

## Résumé des propriétés essentielles

| Propriété | Expression |
|-----------|------------|
| Définition | $T_n(x) = \cos(n \arccos x)$ |
| Récurrence | $T_{n+1}(x) = 2x T_n(x) - T_{n-1}(x)$ |
| Orthogonalité | $\int_{-1}^{1} \frac{T_m T_n}{\sqrt{1-x^2}} dx = \frac{\pi}{2} \delta_{mn} (1 + \delta_{n0})$ |
| Nœuds | $x_k = \cos\left(\frac{2k-1}{2n}\pi\right)$ (racines), $x_k = \cos\left(\frac{k\pi}{n}\right)$ (extrema) |
| Norme infinie | $\max_{x \in [-1,1]} |T_n(x)| = 1$ |
| Dérivée | $T_n'(x) = n U_{n-1}(x)$ |
| Série | $f(x) = \sum_{k=0}^{\infty} a_k T_k(x)$ |
| Coefficients | $a_0 = \frac{1}{\pi} \int_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}} dx$, $a_k = \frac{2}{\pi} \int_{-1}^{1} \frac{f(x) T_k(x)}{\sqrt{1-x^2}} dx$ |

---

## Références

1. Trefethen, L. N. (2000). *Spectral Methods in MATLAB*. SIAM.
2. Boyd, J. P. (2001). *Chebyshev and Fourier Spectral Methods*. Dover Publications.
3. Mason, J. C., & Handscomb, D. C. (2002). *Chebyshev Polynomials*. Chapman & Hall/CRC.
4. Rivlin, T. J. (1990). *Chebyshev Polynomials: From Approximation Theory to Algebra and Number Theory*. Wiley.
```