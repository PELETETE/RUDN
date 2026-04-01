import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve

def chebyshev_points(N, type_points='gauss_lobatto'):
    """Génère les points de collocation de Tchebychev."""
    if type_points == 'gauss_lobatto':
        # Points incluant x = ±1 (Extrema de T_N)
        return np.cos(np.pi * np.arange(N + 1) / N)
    else:
        # Points intérieurs (Racines de T_N)
        return np.cos(np.pi * (2 * np.arange(1, N + 1) - 1) / (2 * N))

def evaluate_chebyshev(n, x):
    """Évalue le n-ième polynôme de Tchebychev T_n(x)."""
    return np.cos(n * np.arccos(np.clip(x, -1, 1)))

def chebyshev_derivative(n, x, order=1):
    """Calcule la dérivée d'ordre 1 ou 2 du n-ième polynôme T_n(x)."""
    x = np.atleast_1d(x)
    if n == 0:
        return np.zeros_like(x)
    if n == 1:
        return np.ones_like(x) if order == 1 else np.zeros_like(x)

    theta = np.arccos(np.clip(x, -1.0, 1.0))
    
    if order == 1:
        # Formule stable pour T_n'(x)
        with np.errstate(divide='ignore', invalid='ignore'):
            res = n * np.sin(n * theta) / np.sin(theta)
        # Gestion des bords x = ±1
        res[np.isnan(res)] = (n**2) * ((-1)**((n+1) * (x[np.isnan(res)] < 0)))
        return res
    
    if order == 2:
        # Formule pour T_n''(x)
        d1 = chebyshev_derivative(n, x, 1)
        t0 = evaluate_chebyshev(n, x)
        with np.errstate(divide='ignore', invalid='ignore'):
            res = (x * d1 - n**2 * t0) / (x**2 - 1)
        # Limite aux bords : T_n''(1) = n^2(n^2-1)/3
        mask = np.isnan(res)
        res[mask] = (n**2 * (n**2 - 1) / 3.0) * ((-1)**(n * (x[mask] < 0)))
        return res

def create_boundary_function(alpha, beta):
    """Crée une fonction affine g(x) telle que g(-1) = alpha et g(1) = beta."""
    return lambda x: 0.5 * (beta * (1 + x) + alpha * (1 - x))

def create_basis_functions(N):
    """Génère les fonctions de base phi_k(x) = T_{k+2}(x) - T_k(x) s'annulant aux bords."""
    basis = []
    basis_deriv = []
    basis_deriv2 = []
    
    for k in range(N):
        def phi(x, k=k):
            return evaluate_chebyshev(k + 2, x) - evaluate_chebyshev(k, x)
        basis.append(phi)
        
        def phi_deriv(x, k=k):
            return chebyshev_derivative(k + 2, x, 1) - chebyshev_derivative(k, x, 1)
        basis_deriv.append(phi_deriv)
        
        def phi_deriv2(x, k=k):
            return chebyshev_derivative(k + 2, x, 2) - chebyshev_derivative(k, x, 2)
        basis_deriv2.append(phi_deriv2)
    
    return basis, basis_deriv, basis_deriv2

def chebyshev_diff_matrix(N):
    """Génère la matrice de dérivation spectrale D."""
    x = chebyshev_points(N, 'gauss_lobatto')
    c = np.ones(N + 1)
    c[0] = c[N] = 2.0
    c = c * (-1.0)**np.arange(N + 1)
    X = np.tile(x, (N + 1, 1))
    dX = X - X.T
    D = (c[:, np.newaxis] / c[np.newaxis, :]) / (dX + np.eye(N + 1))
    D = D - np.diag(np.sum(D, axis=1))
    return D, x



# MÉTHODE DE COLLOCATION 

def solve_bvp_collocation_basis(N, d0, d1, d2, f, alpha, beta, x_plot=None):
    """
    Résout d2(x)*u_xx + d1(x)*u_x + d0(x)*u = f(x) sur [-1,1]
    avec u(-1)=alpha, u(1)=beta
    Méthode de collocation avec fonctions de base
    
    Paramètres:
    -----------
    N : int
        Nombre de fonctions de base
    d0, d1, d2 : callable
        Fonctions coefficients
    f : callable
        Terme source
    alpha, beta : float
        Conditions aux limites u(-1)=alpha, u(1)=beta
    x_plot : array, optional
        Points pour l'évaluation de la solution
    """
    NBASIS = N  # Nombre de fonctions de base
    
    # 1. Points de collocation intérieurs (type Gauss)
    x_colloc = chebyshev_points(NBASIS, 'gauss')
    
    # 2. Fonctions de base et leurs dérivées
    basis, basis_deriv, basis_deriv2 = create_basis_functions(NBASIS)
    
    # 3. Construction du système linéaire H * a = G
    H = np.zeros((NBASIS, NBASIS))
    G = np.zeros(NBASIS)
    
    # Fonction de bord
    boundary_func = create_boundary_function(alpha, beta)
    
    for i in range(NBASIS):
        xi = x_colloc[i]
        
        # Évaluation de la fonction de bord et de sa dérivée
        B = boundary_func(xi)
        Bx = 0.5 * (beta - alpha)  # Dérivée de la fonction affine
        
        # Second membre modifié: G = f - d0*B - d1*Bx
        G[i] = f(xi) - d0(xi) * B - d1(xi) * Bx
        
        # Remplissage de la matrice H
        for j in range(NBASIS):
            phi = basis[j](xi)
            phi_x = basis_deriv[j](xi)
            phi_xx = basis_deriv2[j](xi)
            
            H[i, j] = d2(xi) * phi_xx + d1(xi) * phi_x + d0(xi) * phi
    
    # 4. Résolution du système linéaire
    try:
        a = solve(H, G)
    except np.linalg.LinAlgError:
        print("Erreur : matrice singulière")
        a = np.zeros(NBASIS)
    
    # 5. Reconstruction de la solution
    if x_plot is None:
        x_plot = np.linspace(-1, 1, 101)
    
    u_plot = np.zeros_like(x_plot)
    for i, xi in enumerate(x_plot):
        u_plot[i] = boundary_func(xi)
        for j in range(NBASIS):
            u_plot[i] += a[j] * basis[j](xi)
    
    # Solution aux points de collocation
    u_colloc = np.zeros_like(x_colloc)
    for i, xi in enumerate(x_colloc):
        u_colloc[i] = boundary_func(xi)
        for j in range(NBASIS):
            u_colloc[i] += a[j] * basis[j](xi)
    
    return x_colloc, u_colloc, a, x_plot, u_plot


# MÉTHODE DE COLLOCATION (Version avec matrice de différenciation)


def solve_bvp_collocation_matrix(N, d0, d1, d2, f, alpha, beta):
    """
    Résout d2(x)*u_xx + d1(x)*u_x + d0(x)*u = f(x) sur [-1,1]
    avec u(-1)=alpha, u(1)=beta
    Méthode de collocation avec matrice de différenciation
    """
    # Matrice de différenciation
    D, x = chebyshev_diff_matrix(N)
    D2 = D @ D
    
    # Évaluation des coefficients aux points de collocation
    d0_vals = d0(x)
    d1_vals = d1(x)
    d2_vals = d2(x)
    f_vals = f(x)
    
    # Construction de l'opérateur L = d2*D2 + d1*D + d0*I
    L = np.diag(d2_vals) @ D2 + np.diag(d1_vals) @ D + np.diag(d0_vals)
    
    # Application des conditions aux limites
    # x[0] = 1 (bord droit), x[-1] = -1 (bord gauche)
    L[0, :] = 0
    L[0, 0] = 1
    L[-1, :] = 0
    L[-1, -1] = 1
    
    # Second membre modifié
    F = f_vals.copy()
    F[0] = beta   # u(1) = beta
    F[-1] = alpha # u(-1) = alpha
    
    # Résolution
    try:
        u = solve(L, F)
    except np.linalg.LinAlgError:
        print("Erreur : matrice singulière")
        u = np.zeros(N + 1)
    
    return x, u



# MÉTHODE DE COLLOCATION POUR PROBLÈMES NON-LINÉAIRES

def solve_nonlinear_newton(N, f_nonlinear, alpha, beta, y_init=None, max_iter=50, tol=1e-10):
    """
    Résout y'' = f(x, y, y') avec y(-1)=alpha, y(1)=beta
    par collocation de Tchebychev + méthode de Newton
    
    Paramètres:
    -----------
    N : int
        Nombre de points (incluant les bords)
    f_nonlinear : callable
        Fonction f(x, y, yp) retournant y''
    alpha, beta : float
        Conditions aux limites
    y_init : array, optional
        Initialisation (None = solution linéaire)
    max_iter : int
        Nombre maximum d'itérations Newton
    tol : float
        Tolérance de convergence
    """
    # Matrices de différenciation
    D, x = chebyshev_diff_matrix(N)
    D2 = D @ D
    
    # Initialisation
    if y_init is None:
        # Solution linéaire entre alpha et beta
        y_init = alpha * (1 - x) / 2 + beta * (1 + x) / 2
    
    y = y_init.copy()
    
    for iteration in range(max_iter):
        # Calcul des dérivées
        yp = D @ y
        ypp = D2 @ y
        
        # Résidu: R = ypp - f(x, y, yp)
        R = ypp.copy()
        for i in range(N + 1):
            R[i] -= f_nonlinear(x[i], y[i], yp[i])
        
        # Jacobienne: J = D2 - df/dy - df/dyp * D
        # Approximation par différences finies
        J = D2.copy()
        eps = 1e-8
        
        for i in range(N + 1):
            # Perturbation pour df/dy
            y_pert = y.copy()
            y_pert[i] += eps
            yp_pert = D @ y_pert
            
            f_pert = f_nonlinear(x[i], y_pert[i], yp_pert[i])
            f_current = f_nonlinear(x[i], y[i], yp[i])
            
            df_dy = (f_pert - f_current) / eps
            J[i, i] -= df_dy
        
        # Application des conditions aux limites
        # y(1) = beta (x[0] = 1)
        J[0, :] = 0
        J[0, 0] = 1
        R[0] = y[0] - beta
        
        # y(-1) = alpha (x[-1] = -1)
        J[-1, :] = 0
        J[-1, -1] = 1
        R[-1] = y[-1] - alpha
        
        # Mise à jour Newton
        try:
            delta = solve(J, -R)
        except np.linalg.LinAlgError:
            print(f"Matrice singulière à l'itération {iteration}")
            break
        
        y = y + delta
        
        # Vérification de convergence
        if np.linalg.norm(delta, np.inf) < tol:
            print(f"Convergence Newton en {iteration + 1} itérations")
            break
    
    return x, y


def solve_nonlinear_newton_exact_jacobian(N, f_nonlinear, df_dy, df_dyp,
                                           alpha, beta, y_init=None, max_iter=50, tol=1e-10):
    """
    Version avec jacobienne analytique exacte
    
    Paramètres supplémentaires:
    ---------------------------
    df_dy : callable
        Dérivée partielle ∂f/∂y
    df_dyp : callable
        Dérivée partielle ∂f/∂y'
    """
    D, x = chebyshev_diff_matrix(N)
    D2 = D @ D
    
    if y_init is None:
        y_init = alpha * (1 - x) / 2 + beta * (1 + x) / 2
    
    y = y_init.copy()
    
    for iteration in range(max_iter):
        yp = D @ y
        ypp = D2 @ y
        
        # Résidu
        R = ypp.copy()
        for i in range(N + 1):
            R[i] -= f_nonlinear(x[i], y[i], yp[i])
        
        # Jacobienne analytique
        J = D2.copy()
        for i in range(N + 1):
            dfdy_val = df_dy(x[i], y[i], yp[i])
            dfdyp_val = df_dyp(x[i], y[i], yp[i])
            
            J[i, i] -= dfdy_val
            # Contribution de df/dyp * D
            for j in range(N + 1):
                J[i, j] -= dfdyp_val * D[i, j]
        
        # Conditions aux limites
        J[0, :] = 0
        J[0, 0] = 1
        R[0] = y[0] - beta
        
        J[-1, :] = 0
        J[-1, -1] = 1
        R[-1] = y[-1] - alpha
        
        try:
            delta = solve(J, -R)
        except np.linalg.LinAlgError:
            print(f"Matrice singulière à l'itération {iteration}")
            break
        
        y = y + delta
        
        if np.linalg.norm(delta, np.inf) < tol:
            print(f"Convergence en {iteration + 1} itérations")
            break
    
    return x, y


# EXEMPLE DE TEST


if __name__ == "__main__":
    print("Test des méthodes de collocation de Tchebychev")
    print("="*50)
    
    # Équation: y'' + y = 0, y(-1)=0, y(1)=1
    def d0(x): return 1.0
    def d1(x): return 0.0
    def d2(x): return 1.0
    def f(x): return 0.0
    
    alpha = 0.0
    beta = 1.0
    
    # Solution exacte
    def exact(x):
        return np.sin(x + 1) / np.sin(2)
    
    # Test 1: Méthode avec fonctions de base
    print("\n1. Méthode de collocation avec fonctions de base")
    x_colloc, u_colloc, coeffs, x_plot, u_plot = solve_bvp_collocation_basis(
        20, d0, d1, d2, f, alpha, beta
    )
    error_basis = np.max(np.abs(u_plot - exact(x_plot)))
    print(f"Erreur maximale: {error_basis:.2e}")
    
    # Test 2: Méthode avec matrice de différenciation
    print("\n2. Méthode de collocation avec matrice de différenciation")
    x, u = solve_bvp_collocation_matrix(20, d0, d1, d2, f, alpha, beta)
    u_exact = exact(x)
    error_matrix = np.max(np.abs(u - u_exact))
    print(f"Erreur maximale: {error_matrix:.2e}")
    
    # Affichage des résultats
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.plot(x_plot, exact(x_plot), 'b-', linewidth=2, label='Solution exacte')
    plt.plot(x_colloc, u_colloc, 'ro', markersize=6, label='Points de collocation')
    plt.plot(x_plot, u_plot, 'g--', linewidth=1.5, alpha=0.7, label='Interpolation')
    plt.xlabel('x')
    plt.ylabel('y(x)')
    plt.title('y\'\' + y = 0, y(-1)=0, y(1)=1')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.semilogy(x_plot, np.abs(u_plot - exact(x_plot)), 'r-', linewidth=2)
    plt.xlabel('x')
    plt.ylabel('Erreur')
    plt.title(f'Erreur de la solution (max = {error_basis:.2e})')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # Test non-linéaire: Équation de Bratu
    print("\n" + "="*50)
    print("Test non-linéaire: Équation de Bratu")
    print("y'' + λ exp(y) = 0, y(-1)=0, y(1)=0")
    
    lambda_param = 3.0
    
    def f_bratu(x, y, yp):
        return -lambda_param * np.exp(y)
    
    x, y = solve_nonlinear_newton(30, f_bratu, 0.0, 0.0, max_iter=30)
    
    plt.figure(figsize=(10, 5))
    plt.plot(x, y, 'bo-', linewidth=2, markersize=4)
    plt.xlabel('x')
    plt.ylabel('y(x)')
    plt.title(f'Équation de Bratu: y\'\' + {lambda_param} exp(y) = 0')
    plt.grid(True, alpha=0.3)
    plt.show()
    
    print(f"Valeur maximale de y: {np.max(y):.4f}")