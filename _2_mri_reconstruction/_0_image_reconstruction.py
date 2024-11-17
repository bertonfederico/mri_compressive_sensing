import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage.color import rgb2gray
from skimage.transform import resize
from iterative_algorithms.fft_transforms import ifft2c, Phi
from iterative_algorithms import ISTA_reconstruction, ADMM_TV_reconstruction, ADMM_wavelet_reconstruction
from skimage.metrics import peak_signal_noise_ratio as psnr




##################################################################
######################### Image loading ##########################
##################################################################
# Load and preprocess the input image
file_path = 'image.jpg'  # Path to the input image
image = imread(file_path)

# Convert to grayscale if the image is RGB
if image.ndim == 3:
    image = rgb2gray(image)

# Resize image to desired dimensions (256x256) with anti-aliasing
img_original = resize(image, (256, 256), anti_aliasing=True)



##################################################################
######################### Sampling image #########################
##################################################################
# Gaussian mask creation for k-space
shape = (256, 256)
peak_prob=1
sigma=0.16
center_y, center_x = shape[0] // 2, shape[1] // 2
Y, X = np.ogrid[:shape[0], :shape[1]]
distance_from_center = np.sqrt((X - center_x) ** 2 + (Y - center_y) ** 2)
gaussian_prob = peak_prob * np.exp(- (distance_from_center ** 2) / (2 * (sigma * min(shape)) ** 2))
random_matrix = np.random.rand(*shape)
mask = random_matrix < gaussian_prob  # Sample according to Gaussian probability
gaussian_mask = mask.astype(np.float32)

# Sampling with mask
y_mask_gaussian = Phi(img_original, gaussian_mask)                                              # Undersampled k-space data
img_no_reconstruction_gaussian = ifft2c(y_mask_gaussian)                                        # Direct inverse FFT without reconstruction
mse_ifft_gaussian = np.mean((img_original - np.abs(img_no_reconstruction_gaussian)) ** 2)       # MSE without reconstruction
psnr_ifft_gaussian = psnr(img_original, np.abs(img_no_reconstruction_gaussian))

plt.figure()
plt.suptitle(f"Reconstruction algorithms\nInitial MSE: {mse_ifft_gaussian:.6f}\n Initiali PSNR: {psnr_ifft_gaussian:.6f}", fontsize=25)




##################################################################
#############################  ISTA  #############################
##################################################################
# ISTA reconstruction parameters
lam = 0.01             # Regularization parameter for sparsity
max_iter = 5000        # Maximum number of ISTA iterations
tol = 1e-7             # Convergence tolerance
step_size = 1          # Step size for the ISTA update

# Reconstruct the image using ISTA
img_reconstructed_ista, iterations = ISTA_reconstruction.ISTA_MRI_reconstruction(y_mask_gaussian, gaussian_mask, lam, max_iter, tol, step_size)
mse_ista = np.mean((img_original - np.abs(img_reconstructed_ista)) ** 2)
psnr_insta = psnr(img_original, np.abs(img_reconstructed_ista))

# Print MSE and plot images
recon_type = "ISTA with Wavelet"
variable_setting = f"Lambda: {lam}; Iterations: {iterations}; Step-size: {step_size}"
print(f"\n\n{recon_type}\n{variable_setting}\nMSE: {mse_ista:.8f}\nPSNR: {psnr_insta:.8f}")
plt.subplot(1,3,1)
plt.imshow(np.abs(img_reconstructed_ista))
plt.title(f"{recon_type}\nMSE: {mse_ista:.8f}\nPSNR: {psnr_insta:.8f}")
plt.axis('off')




##################################################################
###########################  ADMM Wavelet ########################
##################################################################
# ISTA reconstruction parameters
wavelet='db1'        # Wavelet type
level=4              # Wavelet level
threshold=0.1        # Threshold
max_iter = 5000      # Maximum number of ADMM iterations

# Reconstruct the image using ADMM
img_reconstructed_admm = ADMM_wavelet_reconstruction.admm_mri_wavelet_reconstruction(y_mask_gaussian, gaussian_mask, num_iters=max_iter)
mse_admm = np.mean((img_original - np.abs(img_reconstructed_admm)) ** 2)
psnr_admm = psnr(img_original, np.abs(img_reconstructed_admm))

# Print MSE and plot images
recon_type = "ADMM with Wavelet"
variable_setting = f"Wavelet: {wavelet} - level {level}; Threshold: {threshold}; Max iterations: {max_iter};"
print(f"\n\n{recon_type}\n{variable_setting}\nMSE: {mse_admm}\nPSNR: {psnr_admm}")
plt.subplot(1,3,2)
plt.imshow(np.abs(img_reconstructed_admm))
plt.title(f"{recon_type}\nMSE: {mse_admm:.8f}\nPSNR: {psnr_admm:.8f}")
plt.axis('off')





##################################################################
############################  ADMM TV  ###########################
##################################################################
# ADMM reconstruction parameters
lam = 0.1             # Regularization parameter for sparsity
max_iter = 100        # Maximum number of ADMM iterations

# Reconstruct the image using ADMM
img_reconstructed_admm = ADMM_TV_reconstruction.admm_mri_tv_reconstruction(y_mask_gaussian, gaussian_mask, lam, max_iter)
mse_admm = np.mean((img_original - np.abs(img_reconstructed_admm)) ** 2)
psnr_admm = psnr(img_original, np.abs(img_reconstructed_admm))

# Print MSE and plot images
recon_type = "ADMM with TV"
variable_setting = f"Lambda: {lam}; Max iterations: {max_iter};"
print(f"\n\n{recon_type}\n{variable_setting}\nMSE: {mse_admm}\nPSNR: {psnr_admm}")
plt.subplot(1,3,3)
plt.imshow(np.abs(img_reconstructed_admm))
plt.title(f"{recon_type}\nMSE: {mse_admm:.8f}\nPSNR: {psnr_admm:.8f}")
plt.axis('off')
plt.show()
