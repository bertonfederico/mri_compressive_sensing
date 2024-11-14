import numpy as np
import pywt
from domain_transforms.phi_functions import Phi, Phi_T
from domain_transforms.fft_transforms import fft2c, ifft2c

# Soft-thresholding function for sparsity regularization
def soft_thresholding(x, threshold):
    """
    Applies soft thresholding to promote sparsity, commonly used in wavelet-based denoising.

    Parameters:
    - x: numpy array or tuple of arrays, the input data to threshold. If a tuple, each element is processed recursively.
    - threshold: float, the threshold level. Values with absolute magnitude below this will be set to zero.

    Returns:
    - Thresholded result, either as a numpy array or a tuple of arrays, with values shrunk towards zero.
    """
    # If x is a tuple (e.g., wavelet coefficient structure), apply soft thresholding to each element
    if isinstance(x, tuple):
        return tuple(soft_thresholding(c, threshold) for c in x)
    
    # Apply soft thresholding for individual values in the array
    return np.sign(x) * np.maximum(np.abs(x) - threshold, 0)


# Function for MRI reconstruction using the ISTA algorithm
def ISTA_MRI_reconstruction(y, mask, lam=0.01, max_iter=100, tol=1e-5, step_size=1, wavelet='db1', level=4):
    """
    Perform MRI image reconstruction using Iterative Shrinkage-Thresholding Algorithm (ISTA).

    Parameters:
    - y: numpy array, observed undersampled k-space data.
    - mask: numpy array, binary mask used for undersampling.
    - lam: float, regularization parameter for controlling sparsity.
    - max_iter: int, maximum number of iterations to perform.
    - tol: float, tolerance for stopping criterion based on convergence.
    - wavelet: str, type of wavelet to use for sparsity in the transform domain.
    - level: int, number of decomposition levels for the wavelet transform.

    Returns:
    - x: numpy array, reconstructed MRI image.
    """
    
    # Initialize the reconstructed image with the zero-filled back-projection of y
    x_or = Phi_T(y, mask)
    x = np.copy(x_or)
    iterations = 0
    

    for k in range(max_iter):
        # Compute the residual in k-space and perform a gradient step
        r = y - Phi(x, mask)
        x_new = x + step_size * Phi_T(r, mask)

        # Wavelet decomposition of x_new and soft-thresholding of coefficients
        coeffs = pywt.wavedec2(x_new, wavelet, level=level)
        coeffs_thresh = [soft_thresholding(c, lam * step_size) for c in coeffs]
        
        # Reconstruct the image from the thresholded wavelet coefficients
        x_new = pywt.waverec2(coeffs_thresh, wavelet)


        # Check for convergence by comparing the norm of the update difference to the tolerance
        if np.linalg.norm(x_new - x) < tol:
            break
        
        # Update the estimate for the next iteration
        x = x_new
        iterations += 1

    return x, iterations
