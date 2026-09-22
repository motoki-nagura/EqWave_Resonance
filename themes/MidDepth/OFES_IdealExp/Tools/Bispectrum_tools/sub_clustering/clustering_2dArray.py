import numpy as np
import matplotlib.pyplot as plt

from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.cluster import DBSCAN


def validate_and_get_pattern_size(patterns):
    patterns = np.asarray(patterns, dtype=float)

    if patterns.ndim != 3:
        raise ValueError(
            "patterns must be a 3D array-like object with shape "
            "(num_patterns, height, width)."
        )

    if patterns.shape[0] == 0:
        raise ValueError("patterns must contain at least one pattern.")

    return patterns, patterns.shape[1:]


def normalize_pattern(pattern):
    pattern = np.asarray(pattern, dtype=float)

    max_abs = np.max(np.abs(pattern))

    if max_abs == 0:
        return pattern

    return pattern / max_abs


def crop_to_shape(pattern, threshold_ratio):
    mask = np.abs(pattern) >= threshold_ratio

    if not np.any(mask):
        return pattern

    rows = np.any(mask, axis=1)
    cols = np.any(mask, axis=0)

    return pattern[rows][:, cols]


def pad_to_size(pattern, size):
    result = np.zeros(size, dtype=float)

    h, w = pattern.shape
    target_h, target_w = size

    copy_h = min(h, target_h)
    copy_w = min(w, target_w)

    result[:copy_h, :copy_w] = pattern[:copy_h, :copy_w]

    return result


def preprocess_pattern(
    pattern,
    size,
    threshold_ratio,
    crop_and_pad=True,
):
    normalized = normalize_pattern(pattern)

    if not crop_and_pad:
        return normalized

    cropped = crop_to_shape(normalized, threshold_ratio)
    padded = pad_to_size(cropped, size)

    return padded


def cluster_patterns(patterns, eps=0.8, threshold_ratio=0.1, crop_and_pad=True,):
    """
    Cluster similar 2D patterns.

    Input:
        patterns:
            A list-like object or numpy.ndarray with shape:
                (num_patterns, height, width)

            Example:
                patterns.shape == (10, 4, 5)

            Values can be any numeric values.
            Each pattern is normalized independently to the range [-1.0, 1.0].

        eps:
            DBSCAN distance threshold.
            Larger values group less similar shapes together.
            Smaller values require shapes to be more similar.

        threshold_ratio:
            Background threshold after normalization.
            Since each pattern is normalized independently,
            values whose absolute values are smaller than threshold_ratio
            are treated as background.

            Example:
                threshold_ratio=0.1 means values smaller than
                10% of that pattern's maximum absolute value
                are treated as background.

        crop_and_pad:
            If True, each normalized pattern is cropped to remove
            surrounding background and then padded back to the
            original size before clustering.

            If False, clustering is performed using only the
            normalized patterns.

    Output:
        groups:
            defaultdict(list)
            Key:
                Cluster label assigned by DBSCAN.
            Value:
                List of pattern indices that belong to that cluster.

        processed_patterns:
            numpy.ndarray with shape:
                (num_patterns, height, width)

            These are the normalized, cropped, and padded patterns
            used as feature vectors for clustering.

        labels:
            numpy.ndarray with shape:
                (num_patterns,)

            Cluster label for each input pattern.

    Note:
        Rotation and reflection are not normalized.
        Therefore, rotated or mirrored shapes are treated as different shapes.
    """

    patterns, size = validate_and_get_pattern_size(patterns)

    processed_patterns = np.array([
        preprocess_pattern(pattern, size, threshold_ratio, crop_and_pad)
        for pattern in patterns
    ])

    features = processed_patterns.reshape(processed_patterns.shape[0], -1)

    model = DBSCAN(
        eps=eps,
        min_samples=1,
        metric="euclidean",
    )

    labels = model.fit_predict(features)

    groups = defaultdict(list)

    for index, label in enumerate(labels):
        groups[label].append(index)

    return groups, processed_patterns, labels


def save_clusters_to_pdf(patterns, groups, processed_patterns, output_pdf):
    patterns = np.asarray(patterns, dtype=float)

    with PdfPages(output_pdf) as pdf:
        for label, indices in groups.items():
            n_patterns = len(indices)

            fig, axes = plt.subplots(
                2,
                n_patterns,
                figsize=(3 * n_patterns, 5),
            )

            if n_patterns == 1:
                axes = np.asarray(axes).reshape(2, 1)

            fig.suptitle(f"Cluster {label}", fontsize=16)

            for col, index in enumerate(indices):
                pattern = patterns[index]

                # Normalize the original pattern just before plotting
                normalized = normalize_pattern(pattern)

                # Show the normalized pattern
                axes[0, col].imshow(
                    normalized,
                    vmin=-1,
                    vmax=1,
                )
                axes[0, col].set_title(f"Normalized #{index}")
                axes[0, col].axis("off")

                # Show the processed pattern used for clustering
                axes[1, col].imshow(
                    processed_patterns[index],
                    vmin=-1,
                    vmax=1,
                )
                axes[1, col].set_title(f"Processed #{index}")
                axes[1, col].axis("off")

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def main():
    patterns = np.array([
        [[0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.8, 0.0, 0.0, 0.0],
         [0.0, 0.9, 0.7, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0]],

        [[0.7, 0.0, 0.0, 0.0, 0.0],
         [0.9, 0.8, 0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0]],

        [[0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, -0.8, 0.0, 0.0, 0.0],
         [0.0, -0.9, -0.7, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0]],

        [[0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.7, 0.8, 0.9, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0]],

        [[0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.6, 0.7, 0.8, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0]],

        [[0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.8, 0.0, 0.0, 0.0],
         [0.0, 0.7, 0.0, 0.0, 0.0],
         [0.0, 0.6, 0.0, 0.0, 0.0]],

        [[0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.9, 0.0, 0.0, 0.0],
         [0.0, 0.8, 0.0, 0.0, 0.0],
         [0.0, 0.7, 0.0, 0.0, 0.0]],

        [[0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.8, 0.8, 0.0, 0.0],
         [0.0, 0.8, 0.8, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0]],

        [[0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.7, 0.7, 0.0, 0.0],
         [0.0, 0.7, 0.7, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0]],

        [[0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, -0.7, -0.7, 0.0, 0.0],
         [0.0, -0.7, -0.7, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0]],
    ])

    groups, processed_patterns, labels = cluster_patterns(
        patterns,
        eps=0.8,
        threshold_ratio=0.1,
        crop_and_pad=True
    )

    print("Labels:")
    print(labels)

    for label, indices in groups.items():
        print(f"Cluster {label}: {indices}")

    save_clusters_to_pdf(
        patterns=patterns,
        groups=groups,
        processed_patterns=processed_patterns,
        output_pdf="clusters.pdf",
    )


if __name__ == "__main__":
    from rich import traceback
    traceback.install()
    main()
