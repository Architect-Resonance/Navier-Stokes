# Formal Analysis: The Discretization Map $\Phi$

## 1. Definition of the Map $\Phi$
Let $\omega(x, t)$ be a smooth vorticity field in $\mathbb{R}^3$. We define a partition of unity $\{ \psi_i \}_{i=1}^N$ subordinate to the support of the intense vorticity (the "vortex filaments").

The map $\Phi: \omega \to G(V, E, W)$ is defined such that:
- **Vertices $V$**: The supports of the basis functions $\psi_i$.
- **Weights $W$**: $w_{ij} = \int_{\mathbb{R}^3} \int_{\mathbb{R}^3} \omega(x) \psi_i(x) \cdot K(x-y) \cdot \omega(y) \psi_j(y) dx dy$
  where $K(x-y)$ is the Biot-Savart kernel.

## 2. The Spectral Bound (Lemma 1)
Let $L_1$ be the Hodge 1-Laplacian on the graph $G$.
**Statement**: There exists a constant $C > 0$ such that the continuous enstrophy dissipation $D(t) = \nu \int |\nabla \omega|^2$ satisfies:
$$ D(t) \geq \nu \cdot \lambda_{\min}(L_1) \cdot \sum_{i} \int |\omega \psi_i|^2 + \text{Error}(\Phi) $$

## 3. Addressing the Circularity of $R$
The invariant $R = 1.85731$ must not be an input. Instead, it must emerge as the **Infimum** of the ratio of spectral gaps for any graph $G$ that can be realized by the map $\Phi$ under the constraint of Lei et al. (directional spanning of $\omega$).

**Revised Conjecture**: 
$$ \inf_{\Phi(\omega) \in \mathcal{S}} \frac{\lambda_{\min}(L_{red})}{\lambda_{\min}(L_{full})} \geq 1.518 $$
where $\mathcal{S}$ is the set of all physically realizable interaction graphs under the enstrophy balance.

## 4. Logical Deductions
1. If $\omega$ satisfies Lei et al., then $G = \Phi(\omega)$ is sufficiently connected (no isolated components).
2. For any such $G$, the Hodge gap $\lambda_{\min}(L_1)$ is bounded by the graph connectivity.
3. The "Snap" at 1.518 represents the minimum dissipative gain required to counteract the maximum possible vortex stretching in the filament model.
