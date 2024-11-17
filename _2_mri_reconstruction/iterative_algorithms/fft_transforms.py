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


def Phi(x, mask):
    """
    Simulates undersampling in k-space by applying a mask to the Fourier transform of an image.

    Parameters:
    - x: numpy array, the input image in the spatial domain.
    - mask: numpy array, binary mask with the same shape as the Fourier-transformed image.

    Returns:
    - undersampled_kspace: numpy array, the masked k-space data after applying the mask.
    """
    # Apply the Fourier transform to move the image to k-space, then mask it for undersampling
    return fft2c(x) * mask


def Phi_T(y, mask):
    """
    Reconstructs an image from undersampled k-space data by applying an inverse Fourier Transform with masking.

    Parameters:
    - y: numpy array, the undersampled k-space data.
    - mask: numpy array, binary mask with the same shape as y, indicating sampled locations.

    Returns:
    - reconstructed_img: numpy array, the reconstructed image in the spatial domain.
    """
    # Mask the k-space data to enforce undersampling, then apply the inverse Fourier transform
    return ifft2c(y * mask)