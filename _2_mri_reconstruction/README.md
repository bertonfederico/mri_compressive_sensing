# Original image reconstruction Using ISTA and ADMM

In this section, after completing the sparse representation of the image and its sampling in the transform domain, the reconstruction of the original image is analyzed. The goal is to obtain an accurate approximation of the source image using several optimization algorithms: ISTA (Iterative Shrinkage-Thresholding Algorithm) and ADMM (Alternating Direction Method of Multipliers), applying them with regularization techniques based on Total Variation and Wavelet Transform. These algorithms make it possible to exploit the sparse properties of the image in the FFT domain to effectively reconstruct it even from undersampled data.

## ISTA (Iterative Shrinkage-Thresholding Algorithm)

The ISTA algorithm is an iterative technique used to solve optimization problems in which an attempt is made to minimize an objective function consisting of two main terms:

- Fitting error (data fidelity term): represents the discrepancy between the estimate and the observed data.
- Regularization: penalizes nonsparsity solutions, that is, solutions in which many variables are nonzero.

The problem then shows itself in the following way:

$$
\min_x J(x) + \lambda \|x\|_1
$$

where:
- $J(*)$ is the loss function,
- $x$ is the variable to be reconstructed,
- $\lambda$ is the regularization parameter that balances the trade-off between the fit to the data and the sparsity of the solution.

The loss function in this case is

$$
J(x) = \frac{1}{2} \| \text{mask} \cdot (\text{FFT2D}(x) - \text{FFT2D}(x_{\text{init}})) \|_2^2
$$

where:
- `mask` represents the sampling mask of the original image,
- $\text{FFT2D}(x_{\text{init}})$ corresponds to the initial frequency matrix.

In the context of the ISTA algorithm, considering the gradient of the cost function to be minimized, it is given by:

$$
\nabla_x J(x) = \text{IFFT2D}(\text{mask} \cdot (\text{FFT2D}(x) - \text{FFT2D}(x_{\text{init}})))
$$

### Extension of ISTA with Wavelet Penalty

The traditional ISTA uses the penalization $\|x\|_1$, i.e., the $l_1$ norm of the image coefficients. However, as already mentioned, the domain that exhibits the best sparsity for this type of images is Wavelet, and can therefore be used to control the image reconstruction.

Assuming the Wavelet transform with an orthonormal matrix (which can be obtained in Python code, for example, via 4-level Haar Wavelet), the solution can be found through soft-thresholding by minimizing the $l_1$ norm. The objective function of ISTA can then be modified to apply the $l_1$ penalization to the Wavelet coefficients, rather than the image values. This leads to the following problem:

$$
\min_x J(x) + \lambda \|W(x)\|_1
$$

### Iterative Steps of the Algorithm

The algorithm can then be implemented in the following way:

1. **Descending gradient step**: the residual $r = y - \Phi(x)$ is calculated, where $\Phi$ represents the transformation related to the FFT and the undersampling mask. Next, a descending gradient step is performed:

$$
x_{\text{temp}}^{(k+1)} = x^k - \alpha \nabla J(x^k) = x^k + \alpha \Phi_T(r)
$$

2. **Soft-thresholding step on Wavelet coefficients**: after the gradient step, the code performs a Wavelet decomposition on the updated $x_{\text{temp}}^{(k+1)}$ image; this transformation decomposes the image into a series of coefficients at different levels of resolution. Next, soft-thresholding is applied to each coefficient: at this point, small coefficients are reduced to zero, leading to a more sparse solution.

$$
z^{(k+1)} = \text{wavelet}(x_{\text{temp}}^{(k+1)})
$$

$$
S_\lambda(z^{(k+1)}) = \text{sign}(z^{(k+1)}) \cdot \max(|z^{(k+1)}| - \lambda, 0)
$$

3. **Image reconstruction step**: once the wavelet coefficients are subject to soft-thresholding, the image is reconstructed through the inverse of the wavelet transform:

$$
x^{(k+1)} = \text{wavelet}^T(S_\lambda(z^{(k+1)}))
$$


4. **Convergence check**: finally, the code checks convergence by comparing the difference between the current solution and the previous solution. If the difference is less than a tolerance "tol", or if the number of iterations reaches a threshold, the algorithm stops.

## ADMM with Total Variation

The Alternating Direction Method of Multipliers (ADMM) is an optimization technique particularly suited for problems that combine multiple cost terms, each requiring a different form of regularization. This approach is well-suited for the reconstruction of undersampled MRI images, where the goal is to simultaneously preserve the certain data and obtain a reconstructed image.

Considering the problem of reconstructing an image $x$ from incomplete k-space data $k$, the objective is to find an image that:

- Respects the measured k-space data (i.e., the sampled frequencies),
- Reduces noise and artifacts that may result from undersampling through regularization.

This methodology applies a thresholding to the Wavelet coefficients of the image. In primal form, the problem is formulated as:

$$
\min_{(x,z)} \frac{1}{2} \| \text{mask} \cdot (\text{FFT2D}(x) - \text{FFT2D}(x_{\text{init}})) \|_2^2 + \lambda \| W(x) \|_1
$$

In augmented Lagrangian form, it is expressed as:

$$
L_\rho(x, z, u) = \frac{1}{2} \| \text{mask} \cdot (\text{FFT2D}(x) - \text{FFT2D}(x_{\text{init}})) \|_2^2 + \lambda \| W(z) \|_1 + \frac{\rho}{2} \| z - x - u \|_2^2
$$

and consequently, via ADMM, it can be solved in three iterative steps:

#### Step 1: $x$ update

In the first step, the resolution of

$$
x^{(k+1)} = \arg \min_x \left( \frac{1}{2} \| \text{mask} \cdot (\text{FFT2D}(x) - \text{FFT2D}(x_{\text{init}})) \|_2^2 + \frac{\rho}{2} \| z^k - x - u^k \|_2^2 \right)
$$

It can be solved by a projection operator that:
- Updates $x$ to minimize $\frac{\rho}{2} \| z^k - x - u^k \|_2^2$, obtaining the solution $x = z^k - u^k$,
- Transforms the matrix values $x$ into the frequency domain, and places the values inside the mask exactly equal to the initial values:

$$
x^{(k+1)} = \text{FFT2D}^{-1} \left( \text{mask} \cdot \text{FFT2D}(x_{\text{init}}) + (1 - \text{mask}) \cdot \text{FFT2D}(z^k - u^k) \right)
$$

#### Step 2: $z$ update

In the second step, we proceed with the update of $z$ in the wavelet domain:

$$
z^{(k+1)} = \arg \min_z \left( \lambda \| W(z) \|_1 + \frac{\rho}{2} \| z - x^{(k+1)} - u^k \|_2^2 \right)
$$

It is possible to calculate the solution by soft-thresholding, having to minimize the norm $l_1$:

$$
z_{\text{wavelet}}^{(k+1)} = \text{softThresholding}_{\frac{\lambda}{\rho}} \left( W(x^{(k+1)} - u^k) \right)
$$

$$
z^{(k+1)} = W^{-1} \left( z_{\text{wavelet}}^{(k+1)} \right) = W^T \left( z_{\text{wavelet}}^{(k+1)} \right)
$$

#### Step 3: $u$ update

Lagrange's multiplier is finally corrected:

$$
u^{(k+1)} = u^k - z^{(k+1)} + x^{(k+1)}
$$

## ADMM with Total Variation

The ADMM algorithm with Total Variation regularization differs from the previous one in the regularization step. Whereas with Wavelet we penalize high-frequency detail coefficients, favoring a sparse representation of structures at multiple scales, TV penalizes the variation of adjacent pixel values, favoring uniform shading.

TV regularization thus minimizes the energy associated with rapid changes in an image. This energy, for an image $x(h, v)$, is defined as:

$$
E^{\text{TV}}(x(h, v)) = \int \sqrt{\left( \frac{\partial x}{\partial h} \right)^2 + \left( \frac{\partial x}{\partial v} \right)^2} dh \ dv
$$

The goal is to find an image $x^*$ that minimizes $E_{\text{TV}}$ but is also faithful to the original image. This involves solving the optimization problem:

$$
x^* = \arg \min_x \left( \frac{1}{2} \| \text{mask} \cdot (\text{FFT2D}(x) - \text{FFT2D}(x_{\text{init}})) \|_2^2 + E^{\text{TV}}(x) \right)
$$

Since $E^{\text{TV}}$ is nonlinear and complex, the problem is solved numerically by the gradient descent method. The algorithm is then implemented iteratively in this way:
- Calculation of the gradients in the vertical and horizontal directions, using finite differences as an estimate.
- Combining the horizontal and vertical variations, adding a stability term $\epsilon$ to make up for the case where the previously calculated differences are zero.
- Calculating the divergence, normalizing the gradients, and calculating the sum of the contributions in the two directions.
- At each iteration, updating image pixels by moving in the direction opposite to the gradient.

In primal-dual form, the problem can be expressed as:

$$
\min_{(x, z)} \frac{1}{2} \left\| \text{mask} \cdot \left( \text{FFT2D}(x) - \text{FFT2D}(x_{\text{init}}) \right) \right\|_2^2 + E^{\text{TV}}(z) \quad \text{s.t.} \quad x = z
$$

Which is solved in an iterative mode via ADMM:
- **Update of $x$**: maximizes the fidelity of the portion of the data that has been kept intact by the mask in the frequency domain (in the same manner as used in the wavelet version):

$$
x^{(k+1)} = \text{FFT2D}^{-1} \left( \text{mask} \cdot \text{FFT2D}(x_{\text{init}}) + (1 - \text{mask}) \cdot \text{FFT2D}(z^k - u^k) \right)
$$

- **Update of $z$**: minimizes the total variance term by performing a TV regularization via a proximal operator, as before:

$$
z^{(k+1)} = \text{proximal}_{\frac{\lambda}{\rho}} \left( x^{(k+1)} + u^k \right)
$$

- **Update of the dual variable $u$**:

$$
u^{(k+1)} = u^k - z^{(k+1)} + x^{(k+1)}
$$

