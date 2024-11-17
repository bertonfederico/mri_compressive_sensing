import numpy as np
import pywt
from .fft_transforms import ifft2c, fft2c

def wavelet_denoise(image, wavelet='db1', level=4, threshold=0.1):
    """
    Applies soft-thresholding denoising to wavelet coefficients.

    Parameters:
        image (numpy.ndarray): Input image to be denoised.
        wavelet (str): Type of wavelet to use for decomposition, e.g., 'db1' (Daubechies).
        level (int): Number of decomposition levels.
        threshold (float): Threshold value for soft-thresholding.

    Returns:
        numpy.ndarray: Denoised image.
    """
    # Decompose the image into wavelet coefficients
    coeffs = pywt.wavedec2(image, wavelet=wavelet, level=level)
    
    # Apply soft-thresholding to the detail coefficients
    coeffs_thresh = [coeffs[0]]  # Keep approximation coefficients unchanged
    for detail_level in coeffs[1:]:
        # Apply soft-thresholding to each tuple of detail coefficients
        coeffs_thresh.append(tuple(pywt.threshold(c, threshold, mode='soft') for c in detail_level))
    
    # Reconstruct the image from the thresholded coefficients
    return pywt.waverec2(coeffs_thresh, wavelet)


def admm_mri_wavelet_reconstruction(kspace_sub, mask, num_iters=100, wavelet='db1', level=4, threshold=0.1):
    """
    Reconstructs an MRI image from undersampled k-space data using ADMM with wavelet-based denoising.

    Parameters:
        kspace_sub (numpy.ndarray): Undersampled k-space data.
        mask (numpy.ndarray): Sampling mask for k-space (1 for sampled, 0 for unsampled).
        num_iters (int): Number of ADMM iterations.
        wavelet (str): Type of wavelet for wavelet denoising.
        level (int): Number of decomposition levels for wavelet transform.
        threshold (float): Threshold for wavelet coefficient denoising.

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
        # Apply wavelet denoising to the real part of (img + u) to encourage sparsity in the wavelet domain
        z = wavelet_denoise(np.real(img + u), wavelet=wavelet, level=level, threshold=threshold)

        # Convert z back to complex type for compatibility with further calculations
        z = z.astype(np.complex128)

        # 3. Dual variable update:
        # Adjust the dual variable u based on the difference between img and z
        # This step enforces consistency between the current estimate and regularized solution
        u += img - z

    # Return the final reconstructed image in the spatial domain
    return img
