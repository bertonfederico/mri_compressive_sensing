import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage.color import rgb2gray
from skimage.transform import resize
from scipy.fftpack import fft2, ifft2, fftshift, ifftshift
from skimage.metrics import mean_squared_error as mse
from skimage.metrics import peak_signal_noise_ratio as psnr





##################################################################
######################### Image loading ##########################
##################################################################
# Load and preprocess the input image
file_path = 'image.jpg'  # Path to the input image
image = imread(file_path)

# Convert the image to grayscale if it is in RGB format
if image.ndim == 3:
    image = rgb2gray(image)

# Resize the image to 256x256
image = resize(image, (256, 256), anti_aliasing=True)

# Perform 2D FFT on the image and shift the zero-frequency component to the center
fft_image = fftshift(fft2(image))





##################################################################
####################### Reconstruct image ########################
##################################################################
def reconstruct_image(sampled_fft):
    """
    Reconstruct the image from its FFT representation.
    
    Args:
        sampled_fft (ndarray): The FFT of the image after applying sampling.
        
    Returns:
        ndarray: The reconstructed image.
    """
    return np.abs(ifft2(ifftshift(sampled_fft)))





##################################################################
####################### Sampling methods #########################
##################################################################
# Percentage of frequency components to retain during sampling
percentage = 10


def gaussian_sampling(fft_image, percentage=50, decay=30):
    """
    Perform Gaussian-based sampling on the FFT image, retaining a percentage of samples
    according to a Gaussian distribution centered at the image's center.
    
    Args:
        fft_image (ndarray): The 2D FFT of the image.
        percentage (int): The percentage of samples to retain.
        decay (int): The standard deviation (decay) of the Gaussian distribution.
        
    Returns:
        sampled_fft (ndarray): The FFT with sampled components.
        mask (ndarray): The mask indicating selected frequency components.
    """
    h, w = fft_image.shape
    y, x = np.indices((h, w))
    center = np.array([h // 2, w // 2])
    
    # Create a 2D Gaussian distribution centered at the middle of the FFT image
    gaussian_prob = np.exp(-((x - center[1]) ** 2 + (y - center[0]) ** 2) / (2 * decay ** 2))
    gaussian_prob /= gaussian_prob.sum()  # Normalize to sum to 1
    
    flattened_prob = gaussian_prob.flatten()
    
    # Calculate the number of samples to keep based on the desired percentage
    total_values = h * w
    num_to_keep = int(total_values * (percentage / 100))
    
    # Create a mask by sampling indices based on Gaussian probabilities
    sampled_indices = np.random.choice(np.arange(total_values), size=num_to_keep, replace=False, p=flattened_prob)
    mask = np.zeros(total_values, dtype=bool)
    mask[sampled_indices] = True
    mask = mask.reshape(h, w)  # Reshape mask to 2D
    
    return fft_image * mask, mask


def random_sampling(fft_image, percentage=50):
    """
    Perform random sampling on the FFT image, retaining a percentage of samples randomly.
    
    Args:
        fft_image (ndarray): The 2D FFT of the image.
        percentage (int): The percentage of samples to retain.
        
    Returns:
        sampled_fft (ndarray): The FFT with randomly selected components.
        mask (ndarray): The mask indicating selected frequency components.
    """
    h, w = fft_image.shape
    total_values = h * w
    num_to_keep = int(total_values * (percentage / 100))
    
    # Generate a random mask with the specified number of selected components
    mask = np.zeros(total_values, dtype=bool)
    selected_indices = np.random.choice(total_values, num_to_keep, replace=False)
    mask[selected_indices] = True
    mask = mask.reshape(h, w)
    
    return fft_image * mask, mask


def threshold_values(fft_image, percentage=50):
    """
    Perform threshold-based sampling on the FFT image, retaining the largest magnitude components.
    
    Args:
        fft_image (ndarray): The 2D FFT of the image.
        percentage (int): The percentage of components with the highest magnitude to retain.
        
    Returns:
        sampled_fft (ndarray): The FFT with selected frequency components based on thresholding.
        mask (ndarray): The mask indicating the selected frequency components.
    """
    magnitude = np.abs(fft_image)  # Get the magnitude of the FFT coefficients
    total_values = magnitude.size
    num_to_keep = int(total_values * (percentage / 100))
    
    # Sort the indices of the FFT coefficients by magnitude (descending order)
    sorted_indices = np.argsort(magnitude.flatten())[::-1]
    threshold_index = sorted_indices[num_to_keep]
    threshold_value = magnitude.flatten()[threshold_index]
    
    # Create a mask to keep only the coefficients greater than or equal to the threshold value
    mask = magnitude >= threshold_value
    
    return fft_image * mask, mask





##################################################################
####################### Apply and plot ###########################
##################################################################
methods = {
    "Random sampling": random_sampling(fft_image, percentage),
    "Gaussian sampling": gaussian_sampling(fft_image, percentage),
    "Top values sampling": threshold_values(fft_image, percentage),
}

# Plot the original image and its FFT magnitude
fig, axes = plt.subplots(1, 2, figsize=(18, 6))
fig.suptitle("Original data", fontsize=16)

# Plot the original image
axes[0].imshow(image, cmap='gray')
axes[0].set_title("Original image")
axes[0].axis('off')

# Plot the FFT magnitude (log scale for better visibility)
axes[1].imshow(np.log1p(np.abs(fft_image)), cmap='plasma')
axes[1].set_title("Original FFT2D magnitude")
axes[1].axis('off')

plt.tight_layout()

# Create subplots to compare the three sampling methods (Random, Gaussian, Top Values)
fig, axes = plt.subplots(3, 3, figsize=(18, 12), gridspec_kw={'height_ratios': [1, 1, 0.3], 'hspace': 0.2})

# Set a title for the entire figure
fig.suptitle("Comparison of different sampling methods (10.0%)", fontsize=16)

# For each sampling method, plot the comparison (mask on top, image sampled below)
for i, (method, (sampled_fft, mask)) in enumerate(methods.items()):
    reconstructed = reconstruct_image(sampled_fft)
    
    # Plot the mask (top part of the window)
    axes[0, i].imshow(mask, cmap='gray')
    axes[0, i].set_title(f"{method} mask")
    axes[0, i].axis('off')
    
    # Plot the reconstructed image (bottom part of the window)
    axes[1, i].imshow(reconstructed, cmap='gray')
    axes[1, i].set_title(f"{method} reconstructed")
    axes[1, i].axis('off')
    
    # Calculate errors (MSE and PSNR)
    error_mse = mse(image, reconstructed)
    error_psnr = psnr(image, reconstructed)
    
    # Plot MSE and PSNR in the third row
    axes[2, i].axis('off')
    axes[2, i].text(0.5, 0.5, f'MSE: {error_mse:.4f}\nPSNR: {error_psnr:.2f} dB', 
                    horizontalalignment='center', verticalalignment='center', fontsize=12)

plt.show()
