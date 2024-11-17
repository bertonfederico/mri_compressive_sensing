# Compression of MRI Images

In the field of medical imaging, particularly with Magnetic Resonance Imaging (MRI), the need for image compression arises from the large amount of data they generate. Indeed, MRI images require a lot of storage space, which can be costly and complicated to manage, especially when dealing with high-resolution scans or multiple sections (slices). Compression techniques therefore become critical to reduce the cost of data storage and transmission, while maintaining sufficient quality for diagnostic purposes.

MRI images are often acquired in grayscale and have some inherent symmetry. This symmetry can be exploited for compression. In particular, an MRI image, representing the spatial domain, can exhibit symmetric properties that allow for data reduction. For example, when an MRI image is grayscale, the left half of the image is often a reflection of the right half. This symmetry suggests that, in some cases, data could be reduced by focusing only on one side of the image. However, symmetry in MRI images is not always perfect, so other compression modes are explored to achieve optimal data size reduction.

## Compression techniques on the MRI k-space

The k-space is the raw data acquired during the MRI scan, which contains all the frequency components needed to reconstruct the image. The idea behind compression in the k-space is that not all frequencies contribute equally to the final image. Some frequencies, particularly those in the center of the k-space, represent the low-frequency components, which generally contain the most significant image information. High-frequency components, located at the edges of the k-space, often represent finer details, which can be more easily removed without significantly compromising image quality.

To explore this idea, we can apply different sampling techniques in the k-space after performing a 2D FFT transform of the MRI image. Three different sampling modes are considered for compression:

1. **Gaussian sampling**: In this mode, k-space values are sampled following a Gaussian distribution, concentrating the sampling more toward the center of the k-space.
   
2. **Random sampling**: In this mode, k-space values are chosen randomly, without any priority, to create a sample selection that does not follow a specific pattern.

3. **Threshold sampling**: here, only the coefficients in the k-space that exceed a certain amplitude threshold are selected, excluding the lowest values, which generally contain less meaningful information for the image.

To compare these modes, the same percentage (10%) of the highest values is taken and the others are set to zero. This approach allows observing the effect of each sampling type on the final image quality and data compression. The results show that the optimal mode is the one that selects the highest coefficients, since it retains the most relevant information for image reconstruction. However, the Gaussian mode does not differ much from it, since the FFT2D in the k-space tends to concentrate mainly around the center, which makes Gaussian sampling similar to sampling the highest values.

![Figure_2](https://github.com/user-attachments/assets/b2363f85-79c4-40f7-99d4-8238c297bbf0)


## Ways of sharing data after FFT2D sampling

Once k-space sampling has been performed, it is important to consider how the data can be shared effectively. The FFT2D of the k-space of an MRI image actually represents the raw form of the data immediately obtained from the MRI. Therefore, it is possible to share the k-space of the MRI directly, but after sampling, it becomes necessary to know the modes by which the data were sampled in order to reconstruct them correctly.

- **Gaussian and random modes**: for these two modes, a fixed sampling mask could be used. In other words, the same sampling mask could be used each time, which would eliminate the need to transmit additional information about how the data were sampled. This solution would simplify data management after sampling.

- **Mode of selecting the most significant values**: In the case of sampling based on selecting the highest values, it is necessary to somehow share additional information regarding the location of the sampled values in the k-space. Since the values are selected non-uniformly, it is necessary to know exactly which coefficients were chosen in order to reconstruct the image correctly.

To address this need, two ways to represent the sampled data were tested:

1. **COO (Coordinate Format)**: this mode represents the data as a sparse matrix, storing only nonzero values and their positions. The index of rows and columns are stored separately, greatly reducing the amount of data to be transmitted compared to a full matrix.

2. **RLE (Run-Length Encoding)**: RLE is a compression technique that represents consecutive sequences of equal values (runs) as a pair of value and frequency (i.e., the value itself and the number of times it repeats). In the context of k-space, this mode can be useful for compressing areas where many coefficients are zero or have similar values, reducing the amount of data to be transmitted.

<table align="center">
  <tr>
    <th style="width: 300px;">Metric</th>
    <th style="width: 200px;">Value</th>
  </tr>
  <tr>
    <td>Original size of FFT2D (in bytes)</td>
    <td>1,048,576</td>
  </tr>
  <tr>
    <td>Size of COO sparse matrix (in bytes)</td>
    <td>157,320</td>
  </tr>
  <tr>
    <td>Size of RLE compressed data (in bytes)</td>
    <td>263,580</td>
  </tr>
  <tr>
    <td>Memory saved percentage (COO)</td>
    <td>85.00%%</td>
  </tr>
  <tr>
    <td>Memory saved percentage (RLE)</td>
    <td>74.86%</td>
  </tr>
</table>

The result is that by applying a sampling of 10 percent of the maximum values in the k-space, only a very small fraction of the original data can be transmitted, i.e., about 15 percent in the case of COO encoding and about 25 percent in the case of RLE encoding.