import numpy as np

# Function to create a Gaussian-based sampling mask in k-space
def generate_gaussian_mask(shape, peak_prob=0.8, sigma=0.2):
    """
    Generates a non-uniform sampling mask based on a Gaussian probability distribution centered in k-space.

    Parameters:
    - shape: tuple of ints (height, width) specifying the mask dimensions.
    - peak_prob: float, the maximum sampling probability at the center of k-space.
    - sigma: float, standard deviation of the Gaussian distribution (0 < sigma <= 1).

    Returns:
    - mask: numpy array, binary mask with non-uniform sampling in k-space.
    """
    # Determine the center coordinates of the mask
    center_y, center_x = shape[0] // 2, shape[1] // 2

    # Generate coordinate grids for Y and X
    Y, X = np.ogrid[:shape[0], :shape[1]]
    
    # Calculate Euclidean distance from each point to the center
    distance_from_center = np.sqrt((X - center_x) ** 2 + (Y - center_y) ** 2)

    # Compute the Gaussian probability map
    gaussian_prob = peak_prob * np.exp(- (distance_from_center ** 2) / (2 * (sigma * min(shape)) ** 2))

    # Generate random values and create the mask based on Gaussian probability
    random_matrix = np.random.rand(*shape)
    mask = random_matrix < gaussian_prob  # Sample according to Gaussian probability

    # Calculate the percentage of elements in the mask
    selected_count = np.sum(mask)
    total_elements = mask.size
    percentage = (selected_count / total_elements) * 100

    return mask.astype(np.float32), percentage
