import numpy as np

# Fourier Transform functions with centered frequency components for MRI applications

def fft2c(img):
    """
    Computes the 2D centered Fourier Transform of an image.

    Parameters:
    - img: numpy array, the input image in the spatial domain.

    Returns:
    - kspace: numpy array, the Fourier-transformed image with frequency components centered.
    """
    # Shift the zero-frequency component to the center, apply 2D Fourier Transform, then shift back
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(img)))

def ifft2c(kspace):
    """
    Computes the 2D inverse centered Fourier Transform to reconstruct an image from k-space.

    Parameters:
    - kspace: numpy array, the input data in k-space (frequency domain).

    Returns:
    - img: numpy array, the reconstructed image in the spatial domain.
    """
    # Shift the zero-frequency component to the center, apply inverse 2D Fourier Transform, then shift back
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(kspace)))