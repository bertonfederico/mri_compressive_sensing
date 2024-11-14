from .fft_transforms import fft2c, ifft2c

# Sampling functions in k-space using a binary mask

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