# MRI signal sparsity

Sparse representation is a fundamental concept in signal processing, particularly useful for compressing and representing signals that contain a large amount of information into a small number of meaningful components. In general, a signal $x$ can be represented in a $\Psi$ basis as a linear combination of a few elements:

$$
x = \Psi s
$$

where:
- $x$ is the original signal.
- $\Psi$ represents the base in which the signal $x$ results scattered.
- $s$ is a vector of sparsity coefficients, i.e., a set of weights that combine the columns of $\Psi$ to reconstruct $x$.


The problem results in optimizing the following function:

$$
\mathbf{s}^* = \arg\min_{\mathbf{s}} \|\mathbf{x} - \Psi \mathbf{s}\|_2^2 \quad \text{s.t.} \quad \|\mathbf{s}\|_0 \leq k,
$$

where $k$ represents the maximum amount of nonzero coefficients you want to use in the sparse form $s$.

The exact resolution of the problem is NP-hard; to overcome computational complexity, a relaxed form of the problem is used, replacing the $\ell_0$ norm (discontinuous and nonconvex) with the convex norm. 

In the case of natural signals, such as MRI, their sparse representation is often well approximated in bases such as Fourier, Wavelet or DCT.

# 2D Fourier transform

In the case of images, the general theory of sparsity can be applied using the Discrete Fourier Transform (DFT), which is an orthonormal basis for representing an image in the frequency domain.

In the case of MRI images, the data acquired are not directly spatial images, but are measurements in what is called the k-space. The k-space is a representation of the spatial frequency of the image, describing how the various frequencies are distributed in the image. To obtain the final image, the data in the k-space must be transformed to the spatial domain by the inverse of the Fourier transform.

The two-dimensional Fourier transform (DFT2D) is a mathematical technique used to transform an image from the spatial domain to the  frequency domain. In mathematical terms, the DFT2D of an image $f[x, y]$ (with $x$ and $y$ representing the spatial coordinates of the image) is given by:

$$
F[u, v] = \sum_{x=0}^{M-1} \sum_{y=0}^{N-1} f[x, y] e^{-j2\pi \left(\frac{ux}{M} + \frac{vy}{N}\right)}
$$

where:
- $F[u, v]$ are the coefficients in the frequency domain,
- $f[x, y]$ is the pixel value in the spatial image,
- $M$ and $N$ are the dimensions of the image.

When an image is transformed into the frequency domain with DFT2D, it is represented as a combination of low frequencies and high frequencies:

- Low frequencies: they are located near the center of the frequency domain and represent global image information, as homogeneous areas;
- High frequencies: these are located further from the center and represent fine details, but also noise.

The 2D DFT is separable along rows and columns, which means it can be calculated as the product of two 1D DFTs. Separable Fourier matrices along rows and columns, therefore, are computed as:

$$
\Psi_{\text{row}}[i, j] = e^{-2\pi i \frac{ij}{n}} \quad \text{e} \quad \Psi_{\text{col}}[i, j] = e^{-2pi i \frac{ij}{m}}
$$

where \(i, j\) are the row and column indices. The image $X$ is then represented as:

$$
X = \Psi_{\text{row}}\ s\  \Psi_{\text{col}}^H
$$

where:
- $\Psi_{\text{row}}$ is the Fourier transform matrix applied to the rows of the image.
- $\Psi_{\text{col}}^H$ is the Hermitian matrix of the Fourier transform applied to the columns.
- $s$ is the vector of the scattered coefficients in the frequency domain.

The signal-to-frequency transformation operation is done through the matrix-coefficient product:

$$
s = \Psi_{\text{row}}^H\ X\ \Psi_{\text{col}}^T
$$

The Python code to perform such a mathematical formulation is given below:

```python
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
```
![Figure_13](https://github.com/user-attachments/assets/c0033707-04f4-4b09-957d-ed848eff6760)


# 2D Wavelet Transform

The Discrete Wavelet Transform (DWT) is a powerful technique for analyzing signals and images, dividing a signal into components at different resolutions. It allows signals to be represented in a sparse manner, concentrating most of the information in low frequencies, and reducing the importance of high frequencies. 

In a one-dimensional context, DWT is used to separate an $x[n]$ signal into two main components:
- approximation(low-frequency), via a low-pass filtering
$$A[n]=∑_k​x[k]h[2n−k]$$
- details (high frequency), via a high-pass filtering
$$D[n]=∑_k​x[k]g[2n−k]$$

The 2D DWT extends the 1D DWT to two-dimensional images (data matrices). In this 2D version, you apply DWT separately to the rows and columns of the image:
1. DWT on rows: each row of the image is subjected to DWT 1D, thus creating a low-frequency approximation matrix and a high-frequency detail matrix for each row
2. DWT on the columns: the 1D DWT is applied separately to the columns of results obtained from the previous step (both approximations and horizontal details). Four sub-bands are then created:
    - LL: vertical approximation in the horizontal approximation matrix (low frequency);
    - LH: vertical approximation in the horizontal detail matrix (high frequency);
    - HL: vertical detail in the horizontal approximation matrix (high frequency);
    - HH: vertical detail in the horizontal detail matrix (high frequency).

By setting an iterative level greater than 1, at each level the LL subband of the previous level is used and is decomposed back into four subbands, preserving the main structure of the signal in LL and obtaining finer and finer details in LH, HL, HH.

In general, DWT leads to a sparse representation because it effectively separates low and high frequencies, allowing only the most significant coefficients to be retained for signal reconstruction.

```python
##################################################################
####################### Wavelet matrix (Ψ) #######################
##################################################################
coeffs = pywt.wavedec2(image, wavelet='db1', level=4)

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
```
![Figure_11](https://github.com/user-attachments/assets/bda5d0b1-90de-4178-8be5-62bb49e10077)


# 2D Discrete Cosine Transform

The Discrete Cosine Transform (DCT) is a transform technique widely used in signal and image processing, especially for compression. The DCT is similar to the Fourier transform, but uses only cosine functions, without the complex Fourier terms. 

The DCT for a sequence of discrete signals $x_0, x_1,..., x_N-1$ is defined as:

$$X_k = \sum_{n=0}^{N-1} x_n \cdot \cos\left(\frac{\pi (2n + 1) k}{2N}\right), \quad k = 0, 1, \dots, N-1
$$

where:
- $x_n$ are the samples of the original signal.
- $X_k$ are the coefficients of the DCT.
- $N$ is the number of samples.
- $k$ is the frequency index in the DCT transform.

In the case of images, the 2D Discrete Cosine Transform is applied separately to the rows and columns of an image.

$$X_{k_1, k_2} = \sum_{n_1=0}^{N_1-1} \sum_{n_2=0}^{N_2-1} x_{n_1, n_2} \cdot \cos\left(\frac{\pi (2n_1 + 1) k_1}{2N_1}\right) \cdot \cos\left(\frac{\pi (2n_2 + 1) k_2}{2N_2}\right)$$

In image processing, DCT (as well as FFT and wavelet) is applied to local blocks (typically of size $8 \cdot 8$) rather than to the entire image. This is because:
- reduces computational complexity: applying DCT on an entire image $M \times N$ requires a computational cost equal to $O(MN \log(MN))$, while instead, by dividing the image into blocks of size $B \times B$, the overall complexity becomes $O(MN \log B)$, with $B \ll \min(M, N)$
- greater sparsity: the application of DCT on local blocks allows a more accurate and sparse representation in each block, thus favoring compressed sensing.


```python

##################################################################
########################## DCT transform #########################
##################################################################
h, w = image.shape
dct_image = np.zeros_like(image)

plt.figure(figsize=(12, 6))
plt.suptitle("DCT sparsity")
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
    plt.title(f"DCT image - block_size {block_size}")
    plt.axis('off')

dct_transform(block_size=4, image_wind=2)
dct_transform(block_size=150, image_wind=3)
dct_transform(block_size=256, image_wind=4)
    
plt.show()
```

![Figure_12](https://github.com/user-attachments/assets/652e4226-5169-408e-bf5b-4d1885490f2e)


## Comparison of sparsity with FFT2D, wavelet and DCT

In the context of the sparse sensinge of image compression, a comparison was made regarding the sparsity obtained through the three techniques used, with the use of different thresholds to establish the negligibility of the values obtained.

As can be seen from the data below, the wavelet has a peculiarity: with a threshold of zero, many of its coefficients are exactly zero, resulting in a relatively high sparsity value even when the threshold is zero.

In contrast, in the FFT and DCT transforms, with threshold equal to zero, no coefficient falls within that range. When the threshold value is increased, however, there is a spike in sparsity: many coefficients are eliminated, leading to a marked reduction in the amount of significant data and increasing sparsity considerably.


| Method    | Threshold   | Below Threshold (%)   | Total Size  
| --------- | ----------- | --------------------- | ------------
| FFT2D     | 0           | 0.00                  | 65536                
| FFT2D     | 0.001       | 67.42                 | 65536                
| Wavelet2D | 0           | 35.46                 | 65536               
| Wavelet2D | 0.001       | 61.69                 | 65536               
| DCT2D     | 0           | 0.00                  | 65536                
| DCT2D     | 0.001       | 71.95                 | 65536                


![Figure_1](https://github.com/user-attachments/assets/fd645493-db7f-4d98-b5e2-0c8b6a5cf3dc)







