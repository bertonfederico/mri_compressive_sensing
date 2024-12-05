# Reconstruction of the original image

In this section, after completing the sparse representation of the image and its sampling in the transform domain, the reconstruction of the original image is analyzed. The goal is to obtain an accurate approximation of the source image using several optimization algorithms: ISTA (Iterative Shrinkage-Thresholding Algorithm) and ADMM (Alternating Direction Method of Multipliers), applying them with regularization techniques based on Total Variation and Wavelet Transform. These algorithms make it possible to exploit the sparse properties of the image in the FFT domain to effectively reconstruct it even from undersampled data.

![Figure_1](https://github.com/user-attachments/assets/f8d1b865-d7f0-4018-8436-036ab6e9f559)

## ISTA (Iterative Shrinkage-Thresholding algorithm)

The ISTA algorithm is an iterative technique used to solve optimization problems in which an attempt is made to minimize an objective function consisting of two main terms:

1. Fitting error (data fidelity term): represents the discrepancy between the estimate and the observed data.
2. Regularization (sparsity penalty): penalizes nonsparsity solutions, that is, solutions in which many variables are nonzero.

The problem then shows itself in the following way:

$$
\min_x\ L(x) + \lambda \| x \|_1
$$

where:
- $L(*)$ is the loss function
- $x$ is the variable to be reconstructed
- $\lambda$ is the regularization parameter that balances the trade-off between the fit to the data and the sparsity of the solution.

The loss function in this case is

$$​L(x)=\frac{1}{2} \|Φ(x)−y\|_2^2$$

where $Φ(x)$ represents the process that maps the image $x$ from the spatial domain to the frequency domain via FFT2D, subsequently applying a mask.

Having to consider in the ISTA algorithm the grdient of the error $r$ to be minimized, we have that

$$∇_x​L(x)=Φ_T(Φ(x)−y)$$

where $Φ_T$ consists of the application of the mask and the subsequent inversion of FFT2D.

### Extension of ISTA with wavelet penalty and $\Phi(x)$ function

Traditional ISTA uses the $\| x \|_1$ penalty, i.e., the $l_1$ norm of the image coefficients. However, as mentioned above the domain with better sparsity is Wavelet.
The use of wavelets makes it possible to represent the image in terms of a set of coefficients that can be sparse. A wavelet transform decomposes an image into a series of coefficients at different resolutions, which describe the structure of the image at different scales.

The ISTA objective function can then be modified to apply the $l_1$ penalty on the wavelet coefficients ($W(x)$), rather than on the image values. This leads to the following problem:

$$\min_x\ \ L(x) + \lambda \| W(x) \|_1$$

### Iterative steps of the algorithm

  1. Descending gradient step: the residual $r = y - \Phi(x)$ is calculated, where $\Phi$ represents the transformation related to the FFT and the undersampling mask. Next, a descending gradient step is performed:

  $$x^{k+1} = x^k - \alpha \nabla L(x^k) = x^k + \alpha \  \Phi_T(r)$$

  2. Soft-thresholding step on Wavelet coefficients: after the gradient step, the code performs a Wavelet decomposition on the updated $x^{k+1}$ image; this transformation decomposes the image into a series of coefficients at different levels of resolution. Next, soft-thresholding is applied to each coefficient: at this point, small coefficients are reduced to zero, leading to a more sparse solution.
 
  $$z^{k+1} = wavelet(x^{k+1})$$
  $$S_\lambda(z^{k+1}) = \text{sign}(z^{k+1}) \cdot \max(|z^{k+1}| - \lambda, 0)$$

  3. Image reconstruction step: Once the wavelet coefficients are subject to soft-thresholding, the image is reconstructed through the inverse of the wavelet transform:
     
  $$x^{k+1} = wavelet_T(S_\lambda(z^{k+1}))$$
  
  3. Convergence check: finally, the code checks convergence by comparing the difference between the current solution and the previous solution. If the difference is less than a tolerance $\text{tol}$, or if the number of iterations reaches a threshold, the algorithm stops. 

### Algorithm code

```python
# Soft-thresholding function for sparsity regularization
def soft_thresholding(x, threshold):
    """
    Applies soft thresholding to promote sparsity, commonly used in Wavelet-based denoising.

    Parameters:
    - x: numpy array or tuple of arrays, the input data to threshold. If a tuple, each element is processed recursively.
    - threshold: float, the threshold level. Values with absolute magnitude below this will be set to zero.

    Returns:
    - Thresholded result, either as a numpy array or a tuple of arrays, with values shrunk towards zero.
    """
    # If x is a tuple (e.g., Wavelet coefficient structure), apply soft thresholding to each element
    if isinstance(x, tuple):
        return tuple(soft_thresholding(c, threshold) for c in x)
    
    # Apply soft thresholding for individual values in the array
    return np.sign(x) * np.maximum(np.abs(x) - threshold, 0)


# Function for MRI reconstruction using the ISTA algorithm
def ISTA_MRI_reconstruction(y, mask, lam, max_iter=1000, tol=1e-5, step_size = 0.001, Wavelet='db1', level=4):
    """
    Perform MRI image reconstruction using Iterative Shrinkage-Thresholding Algorithm (ISTA).

    Parameters:
    - y: numpy array, observed undersampled k-space data.
    - mask: numpy array, binary mask used for undersampling.
    - lam: float, regularization parameter for controlling sparsity.
    - max_iter: int, maximum number of iterations to perform.
    - tol: float, tolerance for stopping criterion based on convergence.
    - Wavelet: str, type of Wavelet to use for sparsity in the transform domain.
    - level: int, number of decomposition levels for the Wavelet transform.

    Returns:
    - x: numpy array, reconstructed MRI image.
    """
    
    # Initialize the reconstructed image with the zero-filled back-projection of y
    x = Phi_T(y, mask)
    

    for k in range(max_iter):
        # Compute the residual in k-space and perform a gradient step
        r = y - Phi(x, mask)
        x_new = x + step_size * Phi_T(r, mask)

        # Wavelet decomposition of x_new and soft-thresholding of coefficients
        coeffs = pywt.wavedec2(x_new, Wavelet, level=level)
        coeffs_thresh = [soft_thresholding(c, lam * step_size) for c in coeffs]
        
        # Reconstruct the image from the thresholded Wavelet coefficients
        x_new = pywt.waverec2(coeffs_thresh, Wavelet)


        # Check for convergence by comparing the norm of the update difference to the tolerance
        if np.linalg.norm(x_new - x) < tol:
            break
        
        # Update the estimate for the next iteration
        x = x_new

    return x
```






## ADMM with TOTAL VARIATION
The Alternating Direction Method of Multipliers (ADMM) is an optimization technique particularly suitable for problems that combine multiple cost terms, each of which requires a different form of regularization. This approach applies well to the reconstruction of undersampled MRI images, where we want to simultaneously preserve the measured data and obtain a “clean” image free of artifacts.
Consider the problem of reconstructing an image x from incomplete k-space data k. The goal is to find an image that:
- respects the measured k-space data (i.e., the observed frequencies)
- reduces noise and artifacts that may result from undersampling through regularization.

In primary form, the problem is posed as.

$$\min_x \frac{1}{2} \| mask \ (FFT2D(x) - k) \|^2_2 + \lambda*TV(x)$$

where:
- $mask$ is the sampling mask, a matrix that retains observed k-space values and ignores unmeasured ones.
- $FFT2D(x)$ is the Fourier transform of the $x$ image, which allows it to be compared with the $k$ data in the k-space domain.
- $TV(x)$ represents the total variation of $x$. Since it cannot be computed in closed form, the minimization of TV is also done iteratively:

  - calculation of the gradient in the horizontal and vertical directions;
  - calculation of the total norm of the gradient and normalization by it of the previously calculated gradients;
  - calculation of both horizontal and vertical divergence, and combination of them to calculate the total divergence for each pixel;
  - updating the image from that divergence, using a desired weight variable;


In primal-dual form, the problem becomes:

$$\min_{x, z} \frac{1}{2} \| mask \ (FFT2D(x) - k) \|^2_2 + \lambda*TV(x)\ \ \ \ \ \ \ s.t\ \ \ x = z$$

which is solved in iterative mode via ADMM:
- update of $x$: minimizes the data fidelity term by keeping $z$ and $u$ fixed.
- update of $z$: minimizes the total variance term by performing a TV regularization.
- updating the dual variable $u$: updates the dual variable to push $x$ and $z$ to coincide, progressively satisfying the $x=z$ constraint.

```python
import numpy as np
from domain_transforms.fft_transforms import fft2c, ifft2c

def _total_variation(image, weight, num_iters=100, tau=0.125):
    """
    Applies denoising using Total Variation (TV) regularization with Gradient Descent.

    Parameters:
        image (numpy.ndarray): Input image to be denoised.
        weight (float): Weight for the TV penalty (higher values increase smoothing).
        num_iters (int): Number of gradient descent iterations.
        tau (float): Gradient descent step size.

    Returns:
        numpy.ndarray: Denoised/reconstructed image.
    """
    # Initialize the working image copy
    img = np.copy(image)
    
    # Iteratively apply gradient descent for TV minimization
    for _ in range(num_iters):
        
        # Calculate gradients in x and y directions (differences between neighboring pixels)
        grad_x = np.roll(img, -1, axis=1) - img  # Gradient along x-axis
        grad_y = np.roll(img, -1, axis=0) - img  # Gradient along y-axis

        # Compute the gradient magnitude for isotropic TV norm
        # Adding a small constant (1e-8) for numerical stability
        grad_norm = np.sqrt(grad_x**2 + grad_y**2 + 1e-8)
        
        # Calculate divergence (TV regularization term)
        # Dividing by grad_norm and computing divergence to achieve smoothness
        div_x = (grad_x / grad_norm) - np.roll(grad_x / grad_norm, 1, axis=1)
        div_y = (grad_y / grad_norm) - np.roll(grad_y / grad_norm, 1, axis=0)
        divergence = div_x + div_y
        
        # Update the image by moving in the direction of the negative gradient
        # The term (tau * weight * divergence) controls the degree of smoothing
        img = img + tau * weight * divergence

    return img


def admm_mri_tv_reconstruction(kspace_sub, mask, lam=0.1, num_iters=50):
    """
    Reconstructs an MRI image from undersampled k-space data using ADMM with Total Variation (TV) regularization.

    Parameters:
        kspace_sub (numpy.ndarray): Undersampled k-space data.
        mask (numpy.ndarray): Sampling mask for k-space (1 for sampled, 0 for unsampled).
        lam (float): Regularization parameter for TV denoising.
        num_iters (int): Number of ADMM iterations.

    Returns:
        numpy.ndarray: Reconstructed image in the spatial domain.
    """
    
    # Initial estimate for the image: inverse Fourier transform of the undersampled k-space
    img = ifft2c(kspace_sub)
    
    # Initialize auxiliary variable z and dual variable u
    z = np.copy(img)
    u = np.zeros_like(img)

    # ADMM Iterative reconstruction loop
    for _ in range(num_iters):
        
        # 1. Data fidelity step:
        # Apply the mask in the Fourier domain to retain sampled values and replace unsampled values
        img = ifft2c(mask * kspace_sub + (1 - mask) * fft2c(z - u))
        
        # 2. Regularization step:
        # Apply TV denoising to the real part of (img + u) to encourage spatial smoothness
        z = _total_variation(np.real(img + u), weight=lam)
        
        # Convert z back to complex type for further calculations
        z = z.astype(np.complex128)
        
        # 3. Dual variable update:
        # Adjust the dual variable u to account for the difference between img and z
        u += img - z

    # Return the final reconstructed image in the spatial domain
    return img
```





## ADMM with Wavelet
The ADMM algorithm with wavelet regularization differs from that with TV regularization in the regularization step. Instead of applying a penalty on the gradients of the image (as in the case of TV), it applies a thresholding to the wavelet coefficients of the image.

In primary form, the problem is posed as.

$$\min_x \frac{1}{2} \| mask*(FFT2D(x) - y) \|^2_2 + \lambda \| W(x) \|_1\ \ \$$

while in primal-dual form

$$\min_{x, z} \frac{1}{2} \| mask*(FFT2D(x) - y) \|^2_2 + \lambda \| z \|_1\ \ \ \ \ \ \ s.t\ \ \ W(x) = z$$

which is solved in iterative mode via ADMM:
- fidelity step, as exhibited in the TV case.
- regularization step: in this case, Wavelet-based denoising is applied. This step exploits the wavelet decomposition of the image, applies thresholding to the detail coefficients, and reconstructs the image using the thresholded coefficients.
- update of the dual variable: as in the case of TV, the dual variable is updated to correct.

The main difference between TV and wavelet regularization lies in the nature of the penalty. While TV penalizes variation in pixel values (favoring uniform shading), wavelet penalizes high-frequency detail coefficients, favoring a sparse representation of structures at multiple scales.

```python
import numpy as np
import pywt
from domain_transforms.fft_transforms import ifft2c, fft2c

def Wavelet_denoise(image, Wavelet='db1', level=4, threshold=0.1):
    """
    Applies soft-thresholding denoising to Wavelet coefficients.

    Parameters:
        image (numpy.ndarray): Input image to be denoised.
        Wavelet (str): Type of Wavelet to use for decomposition, e.g., 'db1' (Daubechies).
        level (int): Number of decomposition levels.
        threshold (float): Threshold value for soft-thresholding.

    Returns:
        numpy.ndarray: Denoised image.
    """
    # Decompose the image into Wavelet coefficients
    coeffs = pywt.wavedec2(image, Wavelet=Wavelet, level=level)
    
    # Apply soft-thresholding to the detail coefficients
    coeffs_thresh = [coeffs[0]]  # Keep approximation coefficients unchanged
    for detail_level in coeffs[1:]:
        # Apply soft-thresholding to each tuple of detail coefficients
        coeffs_thresh.append(tuple(pywt.threshold(c, threshold, mode='soft') for c in detail_level))
    
    # Reconstruct the image from the thresholded coefficients
    return pywt.waverec2(coeffs_thresh, Wavelet)


def admm_mri_Wavelet_reconstruction(kspace_sub, mask, num_iters=100, Wavelet='db1', level=4, threshold=0.1):
    """
    Reconstructs an MRI image from undersampled k-space data using ADMM with Wavelet-based denoising.

    Parameters:
        kspace_sub (numpy.ndarray): Undersampled k-space data.
        mask (numpy.ndarray): Sampling mask for k-space (1 for sampled, 0 for unsampled).
        num_iters (int): Number of ADMM iterations.
        Wavelet (str): Type of Wavelet for Wavelet denoising.
        level (int): Number of decomposition levels for Wavelet transform.
        threshold (float): Threshold for Wavelet coefficient denoising.

    Returns:
        numpy.ndarray: Reconstructed image in the spatial domain.
    """
    
    # Initial estimate for the image: inverse Fourier transform of the undersampled k-space
    img = ifft2c(kspace_sub)
    
    # Initialize auxiliary variable z and dual variable u
    z = np.copy(img)
    u = np.zeros_like(img)

    # ADMM iterative reconstruction loop
    for i in range(num_iters):
        
        # 1. Data fidelity step:
        # Update img by enforcing data consistency in the Fourier domain
        # This step retains sampled k-space values and fills in unsampled areas
        img = ifft2c(mask * kspace_sub + (1 - mask) * fft2c(z - u))

        # 2. Regularization step:
        # Apply Wavelet denoising to the real part of (img + u) to encourage sparsity in the Wavelet domain
        z = Wavelet_denoise(np.real(img + u), Wavelet=Wavelet, level=level, threshold=threshold)

        # Convert z back to complex type for compatibility with further calculations
        z = z.astype(np.complex128)

        # 3. Dual variable update:
        # Adjust the dual variable u based on the difference between img and z
        # This step enforces consistency between the current estimate and regularized solution
        u += img - z

    # Return the final reconstructed image in the spatial domain
    return img
```
