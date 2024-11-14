import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage.color import rgb2gray
from skimage.transform import resize
from masks.random_mask import generate_random_mask
from masks.gaussian_mask import generate_gaussian_mask
from domain_transforms.fft_transforms import ifft2c, fft2c
from domain_transforms.phi_functions import Phi
from iterative_algorithms import ISTA_reconstruction, ADMM_TV_reconstruction, ADMM_wavelet_reconstruction
from mri_images.compression_results import plot_results, print_results




##################################################################
######################### Image loading ##########################
##################################################################
# Load and preprocess the input image
file_path = './mri_images/image.jpg'  # Path to the input image
image = imread(file_path)

# Convert to grayscale if the image is RGB
if image.ndim == 3:
    image = rgb2gray(image)

# Resize image to desired dimensions (256x256) with anti-aliasing
img_original = resize(image, (256, 256), anti_aliasing=True)




##################################################################
######################### Sampling image #########################
##################################################################
# Random mask creation for k-space
random_mask_percentage = 0.5
random_mask = generate_random_mask(img_original.shape, undersampling_factor=random_mask_percentage)
y_random = Phi(img_original, random_mask)                                                  # Undersampled k-space data
img_no_reconstruction_random = ifft2c(y_random)                                            # Direct inverse FFT without reconstruction
mse_ifft_random = np.mean((img_original - np.abs(img_no_reconstruction_random)) ** 2)      # MSE without reconstruction

# Gaussian mask creation for k-space
gaussian_mask, gaussian_mask_percentage = generate_gaussian_mask((256, 256), peak_prob=1, sigma=0.16)
y_gaussian = Phi(img_original, gaussian_mask)                                              # Undersampled k-space data
img_no_reconstruction_gaussian = ifft2c(y_gaussian)                                        # Direct inverse FFT without reconstruction
mse_ifft_gaussian = np.mean((img_original - np.abs(img_no_reconstruction_gaussian)) ** 2)  # MSE without reconstruction



##################################################################
#####################  ISTA - Random mask  #######################
##################################################################
# ISTA reconstruction parameters
lam = 0.01              # Regularization parameter for sparsity
max_iter = 100          # Maximum number of ISTA iterations
tol = 1e-20             # Convergence tolerance
step_size = 1           # Step size for the ISTA update

# Reconstruct the image using ISTA
img_reconstructed_ista, iterations = ISTA_reconstruction.ISTA_MRI_reconstruction(y_random, random_mask, lam, max_iter, tol, step_size)
mse_ista = np.mean((img_original - np.abs(img_reconstructed_ista)) ** 2)

# Print MSE and plot images
recon_type = "ISTA with Wavelet - random mask"
variable_setting = f"Lambda: {lam}; Iterations: {iterations}; Step-size: {step_size}"
print_results(recon_type, variable_setting, mse_ista, mse_ifft_random)
plot_results(img_original, img_no_reconstruction_random, img_reconstructed_ista, mse_ifft_random, mse_ista, random_mask_percentage, random_mask, recon_type)





##################################################################
###################  ISTA -  Gaussian mask #######################
##################################################################
# ISTA reconstruction parameters
lam = 0.01             # Regularization parameter for sparsity
max_iter = 100         # Maximum number of ISTA iterations
tol = 1e-7             # Convergence tolerance
step_size = 1          # Step size for the ISTA update

# Reconstruct the image using ISTA
img_reconstructed_ista, iterations = ISTA_reconstruction.ISTA_MRI_reconstruction(y_gaussian, gaussian_mask, lam, max_iter, tol, step_size)
mse_ista = np.mean((img_original - np.abs(img_reconstructed_ista)) ** 2)

# Print MSE and plot images
recon_type = "ISTA with Wavelet - Gaussian mask"
variable_setting = f"Lambda: {lam}; Iterations: {iterations}; Step-size: {step_size}"
print_results(recon_type, variable_setting, mse_ista, mse_ifft_gaussian)
plot_results(img_original, img_no_reconstruction_gaussian, img_reconstructed_ista, mse_ifft_gaussian, mse_ista, gaussian_mask_percentage, gaussian_mask, recon_type)





##################################################################
#################  ADMM Wavelet - Gaussian mask ##################
##################################################################
# ISTA reconstruction parameters
wavelet='db1'        # Wavelet type
level=4              # Wavelet level
threshold=0.1        # Threshold
max_iter = 500       # Maximum number of ADMM iterations

# Reconstruct the image using ADMM
img_reconstructed_ista = ADMM_wavelet_reconstruction.admm_mri_wavelet_reconstruction(y_gaussian, gaussian_mask, num_iters=max_iter)
mse_ista = np.mean((img_original - np.abs(img_reconstructed_ista)) ** 2)

# Print MSE and plot images
recon_type = "ADMM with Wavelet - Gaussian mask"
variable_setting = f"Wavelet: {wavelet} - level {level}; Threshold: {threshold}; Max iterations: {max_iter};"
print_results(recon_type, variable_setting, mse_ista, mse_ifft_gaussian)
plot_results(img_original, img_no_reconstruction_gaussian, img_reconstructed_ista, mse_ifft_gaussian, mse_ista, gaussian_mask_percentage, gaussian_mask, recon_type)





##################################################################
#####################  ADMM TV - Gaussian mask ###################
##################################################################
# ADMM reconstruction parameters
lam = 0.1             # Regularization parameter for sparsity
max_iter = 100        # Maximum number of ADMM iterations

# Reconstruct the image using ADMM
img_reconstructed_ista = ADMM_TV_reconstruction.admm_mri_tv_reconstruction(y_gaussian, gaussian_mask, lam, max_iter)
mse_ista = np.mean((img_original - np.abs(img_reconstructed_ista)) ** 2)

# Print MSE and plot images
recon_type = "ADMM with TV - Gaussian mask"
variable_setting = f"Lambda: {lam}; Max iterations: {max_iter};"
print_results(recon_type, variable_setting, mse_ista, mse_ifft_gaussian)
plot_results(img_original, img_no_reconstruction_gaussian, img_reconstructed_ista, mse_ifft_gaussian, mse_ista, gaussian_mask_percentage, gaussian_mask, recon_type)


