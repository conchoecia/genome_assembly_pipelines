#!/usr/bin/env python3
"""
Plot Hi-C contact matrix as heatmap.

Reads binned contact matrix from npz file and generates a visualization with
log-scale coloring, chromosome boundaries, and labels.

Usage:
    python plot_hic.py --input matrix.npz --output heatmap.png
"""

import argparse
import sys
import numpy as np
import scipy.sparse
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def load_matrix(matrix_file):
    """
    Load sparse contact matrix and metadata.
    
    Returns:
        tuple: (matrix, chrom_names, bin_size, total_bins)
    """
    print("Loading matrix...", file=sys.stderr)
    
    # Load sparse matrix
    matrix = scipy.sparse.load_npz(matrix_file)
    
    # Load metadata
    metadata_file = matrix_file.replace('.npz', '.metadata.npz')
    metadata = np.load(metadata_file, allow_pickle=True)
    
    chrom_names = metadata['chrom_names']
    bin_size = int(metadata['bin_size'])
    total_bins = int(metadata['total_bins'])
    
    print(f"Matrix shape: {matrix.shape}", file=sys.stderr)
    print(f"Non-zero entries: {matrix.nnz:,}", file=sys.stderr)
    print(f"Bin size: {bin_size:,} bp ({bin_size / 1e6:.2f} Mb)", file=sys.stderr)
    
    return matrix, chrom_names, bin_size, total_bins

def symmetrize_matrix(upper_triangle):
    """
    Mirror upper triangle across diagonal to create symmetric matrix.
    
    Args:
        upper_triangle: sparse matrix in CSR format (upper triangle only)
    
    Returns:
        dense numpy array (full symmetric matrix)
    """
    print("Symmetrizing matrix...", file=sys.stderr)
    
    # Convert to COO for easy manipulation
    coo = upper_triangle.tocoo()
    
    # Create arrays for symmetric matrix
    rows = []
    cols = []
    data = []
    
    for i, j, v in zip(coo.row, coo.col, coo.data):
        rows.append(i)
        cols.append(j)
        data.append(v)
        
        # Add mirror entry if not on diagonal
        if i != j:
            rows.append(j)
            cols.append(i)
            data.append(v)
    
    # Create symmetric sparse matrix
    symmetric = scipy.sparse.coo_matrix(
        (data, (rows, cols)),
        shape=upper_triangle.shape
    )
    
    # Convert to dense for plotting
    print("Converting to dense array...", file=sys.stderr)
    dense = symmetric.toarray()
    
    return dense

def apply_log_transform(matrix):
    """
    Apply log10 transformation to contact counts.
    
    Args:
        matrix: numpy array of contact counts
    
    Returns:
        numpy array of log-transformed values
    """
    print("Applying log10 transformation...", file=sys.stderr)
    
    # Log transform: log10(count + 1)
    log_matrix = np.log10(matrix + 1)
    
    max_val = np.max(log_matrix)
    print(f"Max log10(count): {max_val:.2f}", file=sys.stderr)
    
    return log_matrix

def create_colormap():
    """
    Create custom colormap from white to #B10000.
    
    Returns:
        matplotlib colormap
    """
    # Define colors: white (0) -> dark red (#B10000) (max)
    colors = ['#FFFFFF', '#B10000']
    n_bins = 256
    cmap = LinearSegmentedColormap.from_list('hic', colors, N=n_bins)
    
    return cmap

def get_chromosome_boundaries(chrom_names):
    """
    Calculate chromosome boundary positions for plotting.
    
    Args:
        chrom_names: list of chromosome names (e.g., ['1.0', '1.1', '2.0', '2.1', ...])
    
    Returns:
        tuple: (boundaries, labels, label_positions)
    """
    boundaries = []
    labels = []
    label_positions = []
    
    # Track current position and chromosome groups
    current_pos = 0
    prev_chrom = None
    chrom_start = 0
    
    for i, name in enumerate(chrom_names):
        # Chromosome boundaries occur between different chromosomes
        if prev_chrom is not None:
            curr_chrom = name.rsplit('.', 1)[0]  # Remove .0 or .1
            if curr_chrom != prev_chrom:
                boundaries.append(current_pos)
        
        prev_chrom = name.rsplit('.', 1)[0]
        current_pos += 1
    
    # Calculate label positions (center of each haplotype)
    for i, name in enumerate(chrom_names):
        labels.append(name)
        label_positions.append(i + 0.5)
    
    return boundaries, labels, label_positions

def plot_heatmap(matrix, chrom_names, output_file):
    """
    Generate Hi-C heatmap with chromosome boundaries and labels.
    
    Args:
        matrix: log-transformed contact matrix
        chrom_names: list of chromosome names
        output_file: output PNG file path
    """
    print("Creating figure...", file=sys.stderr)
    
    # Calculate figure size for native resolution (1 pixel per bin)
    dpi = 100
    size_inches = matrix.shape[0] / dpi
    
    fig, ax = plt.subplots(figsize=(size_inches, size_inches), dpi=dpi)
    
    # Create colormap
    cmap = create_colormap()
    
    # Plot heatmap
    print("Plotting heatmap...", file=sys.stderr)
    im = ax.imshow(matrix, cmap=cmap, interpolation='none', aspect='auto', origin='upper')
    
    # Add chromosome boundaries
    boundaries, labels, label_positions = get_chromosome_boundaries(chrom_names)
    
    for boundary in boundaries:
        ax.axhline(boundary, color='black', linewidth=0.5, alpha=0.5)
        ax.axvline(boundary, color='black', linewidth=0.5, alpha=0.5)
    
    # Set axis labels at chromosome centers
    # Only show every other label to avoid crowding
    sparse_labels = [labels[i] if i % 2 == 0 else '' for i in range(len(labels))]
    sparse_positions = label_positions
    
    ax.set_xticks(sparse_positions)
    ax.set_xticklabels(sparse_labels, rotation=90, fontsize=6)
    ax.set_yticks(sparse_positions)
    ax.set_yticklabels(sparse_labels, fontsize=6)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('log10(contact count)', rotation=270, labelpad=15)
    
    # Set title
    ax.set_title('Hi-C Contact Map (Diploid)', fontsize=10, pad=10)
    
    # Tight layout
    plt.tight_layout()
    
    # Save figure
    print(f"Saving to {output_file}...", file=sys.stderr)
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"Saved heatmap to {output_file}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Plot Hi-C contact matrix')
    parser.add_argument('--input', '-i', required=True,
                        help='Input .npz matrix file')
    parser.add_argument('--output', '-o', required=True,
                        help='Output PNG file')
    
    args = parser.parse_args()
    
    # Load matrix
    matrix, chrom_names, bin_size, total_bins = load_matrix(args.input)
    
    # Symmetrize matrix
    symmetric_matrix = symmetrize_matrix(matrix)
    
    # Apply log transformation
    log_matrix = apply_log_transform(symmetric_matrix)
    
    # Plot heatmap
    plot_heatmap(log_matrix, chrom_names, args.output)
    
    print("Done!", file=sys.stderr)

if __name__ == '__main__':
    main()
