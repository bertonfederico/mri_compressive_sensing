import numpy as np

# Function to create a random sampling mask for k-space undersampling
def generate_random_mask(shape, undersampling_factor=0.25):
    """
    Generates a random binary mask for undersampling in k-space.

    Parameters:
    - shape: tuple of ints, the dimensions of the k-space mask (height, width).
    - undersampling_factor: float, fraction of samples to retain (e.g., 0.25 for 25% sampling).

    Returns:
    - mask: numpy array, binary mask with randomly distributed samples according to the specified undersampling factor.
    """
    # Create a random mask where each point has a probability of 'undersampling_factor' to be True
    mask = np.random.rand(*shape) < undersampling_factor

    return mask.astype(np.float32)
