import numpy as np


def find_index_largest_values_NdimArray(arr, ntop, threshold=-np.inf):
    """
    Return the indices of the largest ntop elements in an N-dimensional array.

    Only elements whose values are greater than or equal to threshold are used.

    Parameters
    ----------
    arr : ndarray
        Input array of arbitrary dimension.
    ntop : int
        Number of largest elements to return.
    threshold : float, optional
        Minimum allowed value. Elements smaller than this value are excluded.

    Returns
    -------
    list of tuple
        List of multidimensional indices sorted in descending order
        of the corresponding array values.
    """
    arr = np.asarray(arr)

    if ntop <= 0:
        return []

    flat = arr.ravel()

    # Find flattened indices of elements satisfying the threshold.
    candidates = np.flatnonzero(flat >= threshold)

    if candidates.size == 0:
        return []

    ntop = min(ntop, candidates.size)

    # Select the largest ntop values among the candidates.
    idx = candidates[np.argpartition(flat[candidates], -ntop)[-ntop:]]

    # Sort the selected indices in descending order of the array values.
    idx = idx[np.argsort(flat[idx])[::-1]]

    # Convert flattened indices back to multidimensional indices.
    return [
        tuple(map(int, x))
        for x in zip(*np.unravel_index(idx, arr.shape))
    ]


if __name__ == "__main__":

    np.random.seed(12345)

    # Create a random 4D array.
    arr = np.random.rand(4, 5, 6, 7)

    ntop = 5
    threshold_index = 0.995

    # Find indices of the largest values above the threshold.
    indices = find_index_largest_values_NdimArray(
        arr,
        ntop,
        threshold=threshold_index
    )

    print("Top indices and values:")
    values = []
    for idx in indices:
        value = arr[idx]
        values.append(value)
        print(f"{idx} : {value:.10f}")

    # Reference values obtained by full sorting after threshold filtering.
    reference = np.sort(arr[arr >= threshold_index].ravel())[::-1][:ntop]

    print("\nReference top values:")
    print(reference)

    print("\nValues from returned indices:")
    print(np.array(values))

    if np.allclose(values, reference):
        print("\nTest PASSED.")
    else:
        print("\nTest FAILED.")
