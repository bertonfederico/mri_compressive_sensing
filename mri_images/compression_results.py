import matplotlib.pyplot as plt
import numpy as np
from domain_transforms.fft_transforms import fft2c

def plot_results(img_original, img_no_reconstruction, img_reconstruction, mse_ista, rec_type):
    fig = plt.figure(figsize=(15, 5))
    fig.suptitle(f"{rec_type} reconstruction\nMSE={mse_ista:.10f}", fontsize=20)
    plt.subplot(1, 3, 1)
    plt.imshow(np.abs(img_original), cmap='gray')
    plt.title("Original image")
    plt.axis('off')
    plt.subplot(1, 3, 2)
    plt.imshow(np.abs(img_no_reconstruction), cmap='gray')
    plt.title("Compressed image")
    plt.axis('off')
    plt.subplot(1, 3, 3)
    plt.imshow(np.abs(img_reconstruction), cmap='gray')
    plt.title("Reconstruction")
    plt.axis('off')
    plt.show()


def print_results(recon_type, variable_setting, mse_ista, mse_ifft):
    print("\n\n" + recon_type)
    print(variable_setting)
    print(f"MSE with ISTA reconstruction: {mse_ista}")
    print(f"MSE without reconstruction: {mse_ifft}")


def plot_masks(random_mask, random_mask_percentage, gaussian_mask, gaussian_mask_percentage, img_original):
    fig = plt.figure(figsize=(15, 5))
    fig.suptitle("Compression k-space masks", fontsize=20)
    plt.subplot(2, 3, 2)
    plt.imshow(random_mask)
    plt.title(f"Random mask  -  {random_mask_percentage:.2f} %")
    plt.axis('off')
    plt.axis('off')
    plt.subplot(2, 3, 3)
    plt.imshow(gaussian_mask)
    plt.title(f"Gaussian mask  -  {gaussian_mask_percentage:.2f} %")
    plt.axis('off')
    plt.subplot(2, 3, 4)
    plt.imshow(np.log(np.abs(fft2c(img_original))))
    plt.title(f"Origina FFT2D")
    plt.axis('off')
    plt.subplot(2, 3, 5)
    plt.imshow(np.log(np.abs(fft2c(img_original))) * random_mask)
    plt.title(f"FFT2D with random mask")
    plt.axis('off')
    plt.subplot(2, 3, 6)
    plt.imshow(np.log(np.abs(fft2c(img_original))) * gaussian_mask)
    plt.title(f"FFT2D with Gaussian mask")
    plt.axis('off')
    plt.show()


def plot_compressed_images(real_image, random_image, random_mask_percentage, gaussian_image, gaussian_mask_percentage):
    fig = plt.figure(figsize=(15, 5))
    fig.suptitle("Compressed images", fontsize=20)
    plt.subplot(1, 3, 1)
    plt.imshow(np.abs(real_image), cmap='gray')
    plt.title("Original image")
    plt.axis('off')
    plt.subplot(1, 3, 2)
    plt.imshow(np.abs(random_image), cmap='gray')
    plt.title(f"Random mask  -  {random_mask_percentage:.2f} %")
    plt.axis('off')
    plt.subplot(1, 3, 3)
    plt.imshow(np.abs(gaussian_image), cmap='gray')
    plt.title(f"Gaussian mask  -  {gaussian_mask_percentage:.2f} %")
    plt.axis('off')
    plt.show()
