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
