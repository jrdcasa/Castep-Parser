import numpy as np

def compute_transformation_matrix(box_vectors):
    """
    Compute the transformation matrix for an arbitrary simulation box shape.

    Parameters:
    box_vectors (list): List of box vectors (a, b, c) and angles (alpha, beta, gamma).

    Returns:
    numpy.ndarray: 3x3 transformation matrix.
    """
    a, b, c, alpha, beta, gamma = box_vectors

    # Convert angles to radians
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)

    # Compute box transformation matrix
    cos_alpha = np.cos(alpha_rad)
    cos_beta = np.cos(beta_rad)
    cos_gamma = np.cos(gamma_rad)
    sin_gamma = np.sin(gamma_rad)

    matrix = np.array([
        [a, b*cos_gamma, c*cos_beta],
        [0, b*sin_gamma, c*(cos_alpha - cos_beta*cos_gamma) / sin_gamma],
        [0, 0, c*np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 2*cos_alpha*cos_beta*cos_gamma) / sin_gamma]
    ])

    return matrix

def apply_periodic_bc(positions, transformation_matrix):
    """
    Apply periodic boundary conditions for an arbitrary simulation box shape.

    Parameters:
    positions (numpy.ndarray): N x 3 array of particle positions (x, y, z).
    transformation_matrix (numpy.ndarray): 3x3 transformation matrix for box shape.

    Returns:
    numpy.ndarray: Wrapped particle positions after applying PBC.
    """
    # Transform particle coordinates using the inverse of the transformation matrix
    normalized_positions = np.linalg.solve(transformation_matrix, positions.T).T

    # Apply periodic boundary conditions
    wrapped_positions = normalized_positions % 1.0

    # Transform back to original coordinate space
    wrapped_positions = np.dot(wrapped_positions, transformation_matrix)

    return wrapped_positions


if __name__ == "__main__":
    # Example usage:
    # Define box vectors and angles for a triclinic box (arbitrary shape)
    box_vectors = [10.0, 8.0, 12.0, 90, 90, 90]  # (a, b, c, alpha, beta, gamma)

    # Compute transformation matrix for the triclinic box
    transformation_matrix = compute_transformation_matrix(box_vectors)

    # Generate sample particle positions within the triclinic box
    positions = np.array([[5.0, 4.0, 6.0], [12.0, 8.0, 13.0], [2.0, 10.0, 4.0], [10, 20, 30]])

    # Apply periodic boundary conditions
    wrapped_positions = apply_periodic_bc(positions, transformation_matrix)

    print("Original positions:")
    print(positions)
    print("Wrapped positions (after applying periodic boundary conditions):")
    print(wrapped_positions)