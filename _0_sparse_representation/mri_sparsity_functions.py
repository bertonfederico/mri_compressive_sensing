import numpy as np
import matplotlib.pyplot as plt
import pywt
from skimage.io import imread
from skimage.color import rgb2gray
from skimage.transform import resize
from scipy.fftpack import dct, idct




##################################################################
######################### Image loading ##########################
##################################################################
# Load and preprocess the input image
file_path = '../mri_images/image.jpg'  # Path to the input image
image = imread(file_path)

# Convert to grayscale if the image is RGB
if image.ndim == 3:
    image = rgb2gray(image)

# Resize image to desired dimensions (256x256) with anti-aliasing
img_original = resize(image, (256, 256), anti_aliasing=True)




##################################################################
######################## Fourier transform #######################
##################################################################
n, m = image.shape
omega_n = np.exp(-2j * np.pi / n)
omega_m = np.exp(-2j * np.pi / m)

# Separable Fourier matrix along rows and columns.
Psi_row = np.array([[omega_n ** (i * j) for j in range(n)] for i in range(n)]) / np.sqrt(n)
Psi_col = np.array([[omega_m ** (i * j) for j in range(m)] for i in range(m)]) / np.sqrt(m)

# Sparse matrix: s = Ψ_row^H * x * Ψ_col^T
s = np.dot(Psi_row.conj().T, np.dot(image, Psi_col.T)) 

# Image reconstruction
image_reconstructed = np.real(np.dot(Psi_row, np.dot(s, Psi_col.conj().T)))
plt.figure(figsize=(12, 4))
plt.suptitle("FFT2D sparsity", fontsize=30)
plt.subplot(1, 2, 1)
plt.imshow(image, cmap='gray')
plt.title("Original image")
plt.axis("off")
plt.subplot(1, 2, 2)
plt.imshow(np.log(np.abs(s)))
plt.title(f"log-sparse matrix")
plt.axis("off")
plt.tight_layout()
plt.show()




##################################################################
####################### Wavelet transform ########################
##################################################################
coeffs = pywt.wavedec2(image, wavelet='db1', level=4)

# Creating the figure to display the results
fig, axs = plt.subplots(4, 4, figsize=(12, 12))
plt.suptitle("Wavelet iterations", fontsize=15)
for i in range(1, len(coeffs)):
    # Extract details for each layer
    c = coeffs[i]
    
    # Extract LH, HL, HH
    if isinstance(c, tuple):
        LH, HL, HH = c
    else:
        LH = HL = HH = c

    # LL
    axs[0, i-1].imshow(coeffs[0], cmap='gray')
    axs[0, i-1].set_title(f"LL - Level {i}")
    axs[0, i-1].axis('off')
    # LH
    axs[1, i-1].imshow(LH, cmap='gray')
    axs[1, i-1].set_title(f"LH - Level {i}")
    axs[1, i-1].axis('off')
    # HL
    axs[2, i-1].imshow(HL, cmap='gray')
    axs[2, i-1].set_title(f"HL - Level {i}")
    axs[2, i-1].axis('off')
    # HH
    axs[3, i-1].imshow(HH, cmap='gray')
    axs[3, i-1].set_title(f"HH - Level {i}")
    axs[3, i-1].axis('off')

plt.tight_layout()
plt.show()




##################################################################
########################## DCT transform #########################
##################################################################
h, w = image.shape
dct_image = np.zeros_like(image)

plt.figure(figsize=(12, 6))
plt.suptitle("DCT sparsity", fontsize=30)
plt.subplot(1, 4, 1)
plt.imshow(image, cmap='gray')
plt.title("Original image")
plt.axis('off')

# Block iterations
def dct_transform(block_size, image_wind):
    for i in range(0, h, block_size):
        for j in range(0, w, block_size):
            # Extract the block
            block = image[i:i+block_size, j:j+block_size]
            # Apply DCT to the block
            dct_image[i:i+block_size, j:j+block_size] = dct(dct(block.T, norm='ortho').T, norm='ortho')
    plt.subplot(1, 4, image_wind)
    plt.imshow(np.log(np.abs(dct_image) + 1), cmap='hot')
    sparsity = 1 - np.count_nonzero(dct_image) / dct_image.size
    plt.title(f"DCT image - block_size {block_size}\nSparsity: {sparsity:.2%}")
    plt.axis('off')

dct_transform(block_size=4, image_wind=2)
dct_transform(block_size=128, image_wind=3)
dct_transform(block_size=256, image_wind=4)
plt.show()
