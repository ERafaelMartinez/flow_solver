import matplotlib.pyplot as plt
import random
import numpy as np

def plot_data_samples(inputs, outputs, title="", sample_idx=None):
    """
    Plots the input and output data for a specific sample index.
    
    Args:   
        inputs: Input data (torch.Tensor or np.ndarray)
        outputs: Output data (torch.Tensor or np.ndarray)
        title: Title of the plot
        sample_idx: Index of the sample to plot. If None, a random index is chosen.
    """
    if sample_idx is None:
        sample_idx = random.randint(0, len(inputs) - 1)

    # Helper function to convert to numpy if needed
    def to_numpy(data):
        if hasattr(data, "numpy"):
            return data.numpy()
        return data

    input_sample = to_numpy(inputs[sample_idx])
    output_sample = to_numpy(outputs[sample_idx])

    # Collect fields to plot. List of tuples: (type, data, label)
    plots_to_show = []

    # corresponds to a horizontal velocity field (Input u)
    plots_to_show.append(('scalar', input_sample[0], f"Sample Input u (Idx: {sample_idx})"))

    # corresponds to a vertical velocity field (Input v, if present)
    if input_sample.shape[0] > 1:
        plots_to_show.append(('scalar', input_sample[1], f"Sample Input v (Idx: {sample_idx})"))

    # corresponds to a horizontal velocity field (Output u)
    plots_to_show.append(('scalar', output_sample[0], f"Sample Output u (Idx: {sample_idx})"))

    # corresponds to a vertical velocity field (Output v, if present)
    if output_sample.shape[0] > 1:
        plots_to_show.append(('scalar', output_sample[1], f"Sample Output v (Idx: {sample_idx})"))

    # Add velocity field quiver plots if applicable (Input)
    if input_sample.shape[0] >= 2:
        plots_to_show.append(('vector', input_sample, f"Input Velocity Field (Idx: {sample_idx})"))
    
    # Add velocity field quiver plots if applicable (Output)
    if output_sample.shape[0] >= 2:
        plots_to_show.append(('vector', output_sample, f"Output Velocity Field (Idx: {sample_idx})"))

    # Compute global vmin/vmax for consistent coloring of scalar plots
    scalar_data = [p[1].flatten() for p in plots_to_show if p[0] == 'scalar']
    if scalar_data:
        all_data = np.concatenate(scalar_data)
        vmin = all_data.min()
        vmax = all_data.max()
    else:
        vmin, vmax = 0, 1 # Fallback, though unlikely to be reached

    n_plots = len(plots_to_show)
    # Adjust figsize to accommodate potentially more plots (scalars + vectors)
    # Using a base width per plot to ensure readability
    fig, axes = plt.subplots(1, n_plots, figsize=(5 * n_plots, 5), constrained_layout=True)
    
    # Ensure axes is iterable even if only 1 plot
    if n_plots == 1:
        axes = [axes]

    im = None
    scalar_axes = []
    
    for ax, (ptype, data, label) in zip(axes, plots_to_show):
        if ptype == 'scalar':
            im = ax.imshow(data, origin='lower', cmap="RdYlBu_r", vmin=vmin, vmax=vmax)
            ax.set_title(label)
            scalar_axes.append(ax)
        elif ptype == 'vector':
            # data shape: (C, H, W)
            h, w = data.shape[1], data.shape[2]
            x, y = np.arange(w), np.arange(h)
            X, Y = np.meshgrid(x, y)
            
            u = data[0]
            v = data[1]
            magnitude = np.sqrt(u**2 + v**2)
            
            quiv = ax.quiver(X, Y, u, v, magnitude, cmap='viridis', angles='xy', scale_units='xy',
            scale=magnitude.max()/1.5)
            fig.colorbar(quiv, ax=ax, label='Velocity Magnitude')
            
            ax.set_title(label)
            ax.set_aspect('equal')
            ax.grid(True, linestyle=':', alpha=0.6)

    # Add single horizontal colorbar below the scalar plots
    # Only attach it to the scalar axes
    if scalar_axes and im is not None:
        fig.colorbar(im, ax=scalar_axes, orientation='horizontal', fraction=0.05, aspect=50, label="Value")

    plt.suptitle(title)
    plt.show()
