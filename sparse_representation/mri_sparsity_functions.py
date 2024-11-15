import numpy as np
import matplotlib.pyplot as plt
import pywt
from skimage.io import imread
from skimage.color import rgb2gray
from skimage.transform import resize




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


'''

##################################################################
####################### Fourier matrix (Ψ) #######################
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
plt.figure(figsize=(12, 6))
plt.suptitle("FFT2D sparsity")
plt.subplot(1, 3, 1)
plt.imshow(image, cmap='gray')
plt.title("Original image")
plt.axis("off")
plt.subplot(1, 3, 2)
plt.imshow(np.log(np.abs(s)))
plt.title("log-sparse matrix")
plt.axis("off")
plt.subplot(1, 3, 3)
plt.imshow(image_reconstructed, cmap='gray')
plt.title("Reconstructed image")
plt.axis("off")
plt.tight_layout()
plt.show()

total_coeffs = s.size  # Il numero totale di coefficienti in s

# Calcolare il numero di coefficienti non nulli (sparsità)
epsilon = 0.1  # Soglia per considerare un coefficiente come nullo
non_zero_coeffs = np.sum(np.abs(s) > epsilon)  # Conta i coefficienti non nulli

# Calcolare la sparsità
sparsity = non_zero_coeffs / total_coeffs
print(sparsity)


'''

##################################################################
####################### Wavelet matrix (Ψ) #######################
##################################################################
coeffs = pywt.wavedec2(image, wavelet='db1', level=4)





coeffs2 = pywt.wavedec2(image, 'db1', level=4)

# Estrai il coefficiente LL e i dettagli LH, HL, HH per ciascun livello
LL = coeffs2[0]  # Bassa frequenza, il primo elemento
detail_coeffs = coeffs2[1:]  # Dettagli per ciascun livello

# Calcola il numero totale di coefficienti (somma delle dimensioni)
total_coeffs = sum([c.size for c in [LL] + [item for sublist in detail_coeffs for item in sublist]])

# Calcola il numero di coefficienti nulli (con una soglia)
epsilon = 0.1  # Soglia per considerare un coefficiente nullo
null_coeffs = np.sum(np.abs(LL) < epsilon)  # Coefficienti nulli per LL
for level in detail_coeffs:
    for coeff in level:
        null_coeffs += np.sum(np.abs(coeff) < epsilon)  # Somma dei coefficienti nulli per ogni livello

# Calcola la sparsità
sparsity = null_coeffs / total_coeffs

# Mostra i risultati
print(f"Numero totale di coefficienti: {total_coeffs}")
print(f"Numero di coefficienti nulli: {null_coeffs}")
print(f"Sparsità: {sparsity:.4f}")

# Creating the figure to display the results
fig, axs = plt.subplots(4, 4, figsize=(12, 12))

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

