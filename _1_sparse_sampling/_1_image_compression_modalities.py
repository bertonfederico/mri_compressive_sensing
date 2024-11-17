import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft2, fftshift
from scipy.sparse import coo_matrix
from skimage.io import imread
from skimage.color import rgb2gray
from skimage.transform import resize

# Load and resize the image
file_path = 'image.jpg'
image = imread(file_path)

# Convert the image to grayscale if it is in RGB format
if image.ndim == 3:
    image = rgb2gray(image)

# Resize the image to 256x256
image_resized = resize(image, (256, 256), anti_aliasing=True)

# Perform 2D FFT and apply the shift
fft_image = fft2(image_resized)
fft_image_shifted = fftshift(fft_image)

# Select the top 10% of the highest FFT2D values
magnitude = np.abs(fft_image_shifted)
threshold = np.percentile(magnitude, 90)

# Set all values below the threshold to zero
fft_image_shifted[magnitude < threshold] = 0

# Compress using COO (Coordinate Format)
rows, cols = np.nonzero(fft_image_shifted)
values = fft_image_shifted[rows, cols]

# Create the sparse matrix in COO format
sparse_matrix = coo_matrix((values, (rows, cols)), shape=fft_image_shifted.shape)

# Function for Run-Length Encoding (RLE) compression
def rle_compress(values):
    """Compress the non-zero values using Run-Length Encoding."""
    if len(values) == 0:
        return []
    
    compressed = []
    prev_value = values[0]
    count = 1
    
    # Iterate through the values and compress consecutive identical values
    for value in values[1:]:
        if value == prev_value:
            count += 1
        else:
            compressed.append((prev_value, count))
            prev_value = value
            count = 1
    
    # Append the last run
    compressed.append((prev_value, count))
    
    return compressed

# Compress the non-zero values using RLE
rle_compressed = rle_compress(values)

# Calculate the memory usage of the COO format
coo_size = sparse_matrix.data.nbytes + sparse_matrix.row.nbytes + sparse_matrix.col.nbytes

# Calculate the memory usage of the compressed data using RLE
rle_size = sum(len(str(value)) + len(str(count)) for value, count in rle_compressed)

# Calculate the memory usage of the original FFT2D
original_size = fft_image_shifted.nbytes  # Memory occupied by the original FFT2D

# Output the results
print(f"Original size of FFT2D (in bytes): {original_size}")
print(f"Size of COO sparse matrix (in bytes): {coo_size}")
print(f"Size of RLE compressed data (in bytes): {rle_size}")
print(f"Memory saved percentage (COO): {(1-(original_size - coo_size)/original_size)*100:.2f}%")
print(f"Memory saved percentage (RLE): {(1-(original_size - rle_size)/original_size)*100:.2f}%")
