# Compressive sensing for MRI 🧠🔍

This project aims to study the compression of magnetic resonance imaging (MRI) images, focusing on the mode in which these images are initially acquired, namely in the k-space format. The k-space represents a complex matrix that contains spatial information of the acquired images, but is not directly interpretable as a visible image. This type of data is subsequently transformed, via a Fourier transform, into a spatial image. The main goal of the project is to develop techniques to effectively compress MRI images starting from k-space and trying to maintain good visual quality with reduced data volume. The project has been divided into three main sub-projects: the initial study of medical image sparsity, FFT2 sampling of MRI images, and finally optimization techniques for approximate reconstruction of the original images.

## Study of sparsity in medical images 💡📉
The first sub-project explores possible ways to generate sparsity in medical images. Specifically, several transformations, such as two-dimensional FFT (FFT2D), wavelets, and the Discrete Cosine Transform (DCT), were tested, with the goal of identifying which of these techniques can best represent medical images in sparse form. 
- Link for developed code: [tap here](https://github.com/bertonfederico/mri_compressive_sensing/tree/a7b0b243e7640aa5545b67505fa9105cf5a8a54e/_0_sparse_representation)
- Link for related README.md file: [tap here](https://github.com/bertonfederico/mri_compressive_sensing/blob/a7b0b243e7640aa5545b67505fa9105cf5a8a54e/_0_sparse_representation/README.md)

## FFT2D sampling for MRI images 🎯📊
The second sub-project focuses on FFT2D, as MRI images are already represented in this format. In particular, the sampling methods of this domain are studied to reduce the amount of data needed. Several types of sampling masks were explored: a random mask, a centered Gaussian mask and a mask that selects values in the FFT2D that exceed a certain amplitude threshold. In addition, two encoding modes for transmitting sampled data were explored: COO (Coordinate List) and RLE (Run-Length Encoding).
- Link for developed code: [tap here](https://github.com/bertonfederico/mri_compressive_sensing/blob/a7b0b243e7640aa5545b67505fa9105cf5a8a54e/_1_sparse_sampling)
- Link for related README.md file: [tap here](https://github.com/bertonfederico/mri_compressive_sensing/blob/a7b0b243e7640aa5545b67505fa9105cf5a8a54e/_1_sparse_sampling/README.md)

## Optimization techniques for image reconstruction 🧩🔄
The third sub-project focuses on the optimization techniques used to reconstruct or approximate the original image from the compressed data. Two main algorithms were explored: ISTA (Iterative Shrinkage Thresholding Algorithm) and ADMM (Alternating Direction Method of Multipliers), with the aim of reducing the reconstruction error. In addition, two regularization techniques were used: the Wavelet transform and Total Variation, which aim to preserve the main features of the image during the reconstruction process. 
- Link for developed code: [tap here](https://github.com/bertonfederico/mri_compressive_sensing/blob/8bf0a76861771bff9b10f03cef6cab3ca6a90f95/_2_mri_reconstruction)
- Link for related README.md file: [tap here](https://github.com/bertonfederico/mri_compressive_sensing/blob/8bf0a76861771bff9b10f03cef6cab3ca6a90f95/_2_mri_reconstruction/README.md)
