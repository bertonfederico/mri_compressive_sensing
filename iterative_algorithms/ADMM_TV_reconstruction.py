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
