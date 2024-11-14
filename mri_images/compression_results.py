import matplotlib.pyplot as plt
import numpy as np
from domain_transforms.fft_transforms import fft2c

def plot_results(img_original, img_no_reconstruction, img_reconstruction, mse_ifft, mse_ista, percentage, mask, rec_type):
    fig = plt.figure(figsize=(15, 5))
    fig.suptitle(rec_type)
    plt.subplot(2, 3, 1)                                      # Original image
    plt.imshow(np.abs(img_original), cmap='gray')
    plt.title("Original image")
    plt.axis('off')
    plt.subplot(2, 3, 2)                                      # FFT of original image
    plt.imshow(np.log(np.abs(fft2c(img_original))))
    plt.title("log-FFT of original image")
    plt.axis('off')
    plt.subplot(2, 3, 3)                                      # k-space mask
    plt.imshow(mask, cmap='gray')
    plt.title(f"k-Space mask\nUndersampling factor={percentage:.2f}")
    plt.axis('off')
    plt.subplot(2, 3, 4)                                      # Image without reconstruction
    plt.imshow(np.abs(img_no_reconstruction), cmap='gray')
    plt.title(f"Without reconstruction\nMSE={mse_ifft:.10f}")
    plt.axis('off')
    plt.subplot(2, 3, 5)                                      # ISTA/ADMM reconstruction
    plt.imshow(np.abs(img_reconstruction), cmap='gray')
    plt.title(f"Reconstruction\nMSE={mse_ista:.10f}")
    plt.axis('off')
    plt.show()


def print_results(recon_type, variable_setting, mse_ista, mse_ifft):
    print("\n\n" + recon_type)
    print(variable_setting)
    print(f"MSE with ISTA reconstruction: {mse_ista}")
    print(f"MSE without reconstruction: {mse_ifft}")
