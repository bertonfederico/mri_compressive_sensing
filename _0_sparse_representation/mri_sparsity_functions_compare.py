import numpy as np
import matplotlib.pyplot as plt
import pywt
from skimage.io import imread
from skimage.color import rgb2gray
from skimage.transform import resize
from scipy.fftpack import dct, fft2
import time





##################################################################
######################### Image Loading ##########################
##################################################################
# Load and preprocess the input image
file_path = '../mri_images/image.jpg'  # Path to the input image
image = imread(file_path)

# Convert to grayscale if the image is RGB
if image.ndim == 3:
    image = rgb2gray(image)

# Resize image to desired dimensions (256x256) with anti-aliasing
image = resize(image, (256, 256), anti_aliasing=True)

# Thresholds for sparsity calculation
thresholds = [0, 0.001]





##################################################################
########################### FFT2D ###############################
##################################################################
# Compute FFT2D using scipy fft2
fft_result = fft2(image)

# Normalize absolute values between 0 and 1
fft_abs_normalized = np.abs(fft_result) / np.max(np.abs(fft_result))

# Compute percentage under threshold for normalized absolute values
fft_percent_below_threshold = {}
for t in thresholds:
    fft_percent_below_threshold[t] = (
        np.sum(fft_abs_normalized <= t) / fft_result.size * 100
    )

fft_total_size = fft_result.size





##################################################################
########################### Wavelet ############################
##################################################################
# Compute Wavelet transform using PyWavelets
coeffs = pywt.wavedec2(image, wavelet='db1', level=4)

# Concatenate all wavelet coefficients (LL, LH, HL, HH at all levels)
wavelet_total_matrix = []
for i in range(1, len(coeffs)):
    c = coeffs[i]
    if isinstance(c, tuple):
        wavelet_total_matrix.extend([c[0].flatten(), c[1].flatten(), c[2].flatten()])
    else:
        wavelet_total_matrix.append(c.flatten())
wavelet_total_matrix.append(coeffs[0].flatten())
wavelet_total_matrix = np.concatenate(wavelet_total_matrix)

# Normalize absolute values between 0 and 1
wavelet_abs_normalized = np.abs(wavelet_total_matrix) / np.max(np.abs(wavelet_total_matrix))

# Compute percentage under threshold for normalized absolute values
wavelet_percent_below_threshold = {}
for t in thresholds:
    wavelet_percent_below_threshold[t] = (
        np.sum(wavelet_abs_normalized <= t) / wavelet_total_matrix.size * 100
    )

wavelet_total_size = wavelet_total_matrix.size
print(wavelet_total_size)






##################################################################
########################### DCT ###############################
##################################################################
# Compute DCT
dct_result = dct(dct(image.T, norm='ortho').T, norm='ortho')  # Full-image DCT

# Normalize absolute values between 0 and 1
dct_abs_normalized = np.abs(dct_result) / np.max(np.abs(dct_result))

# Compute percentage under threshold for normalized absolute values
dct_percent_below_threshold = {}
for t in thresholds:
    dct_percent_below_threshold[t] = (
        np.sum(dct_abs_normalized <= t) / dct_result.size * 100
    )

dct_total_size = dct_result.size





##################################################################
######################### Visualization ##########################
##################################################################
# Visualize the results of the three transforms (no threshold applied here)
fig, axs = plt.subplots(1, 3, figsize=(18, 4))
fig.suptitle("Comparison: FFT2D, Wavelet-LL, DCT", fontsize = 20)

# FFT
axs[0].imshow(np.log(np.abs(fft_result) + 1), cmap='hot')
axs[0].set_title("FFT2D")
axs[0].axis('off')

# Wavelet (show only LL component as example)
axs[1].imshow(np.log(np.abs(coeffs[0]) + 1), cmap='hot')
axs[1].set_title("Wavelet (LL Component) with 4 levels")
axs[1].axis('off')

# DCT
axs[2].imshow(np.log(np.abs(dct_result) + 1), cmap='hot')
axs[2].set_title("DCT")
axs[2].axis('off')

plt.tight_layout()
plt.show()

##################################################################
######################### Summary Table ##########################
##################################################################
# Helper function to format time in scientific notation
def format_time_scientific(t):
    return f"{t:.2e}"

# Generate a comparison table for the three transforms
print("\nComparison between the three methods (Normalized Absolute Values Below Threshold):")

header = f"{'Method':<10}{'Threshold':<12}{'Below Threshold (%)':<20}{'Total Size':<12}"
print(header)
print("-" * len(header))

for method, size, percent_below in [
    ("FFT2D", fft_total_size, fft_percent_below_threshold),
    ("Wavelet", wavelet_total_size, wavelet_percent_below_threshold),
    ("DCT", dct_total_size, dct_percent_below_threshold),
]:
    for t in thresholds:
        print(f"{method:<10}{t:<12}{percent_below[t]:<20.2f}{size:<12}")
