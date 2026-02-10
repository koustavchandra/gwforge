import numpy as np

# Exact index-based sampling (discrete, no interpolation error)
class Interped2D_Discrete(Prior):
    """Index-based 2D Prior: sample grid indices exactly, optionally jitter within cells."""
    def __init__(self, name, x_axis, y_axis, pdf_grid, x_name='x', y_name='y', jitter=True):
        super().__init__(name=name, latex_label=f"${name}$", unit=None)
        self.x_axis = np.asarray(x_axis)
        self.y_axis = np.asarray(y_axis)
        self.x_name = x_name
        self.y_name = y_name
        self.jitter = jitter

        # Normalize grid to PMF (discrete probabilities)
        self.pmf = np.asarray(pdf_grid).flatten() / np.sum(pdf_grid)
        self.nx, self.ny = len(x_axis), len(y_axis)
        
        # Cell widths for jittering
        self.dx = x_axis[1] - x_axis[0] if len(x_axis) > 1 else 1.0
        self.dy = y_axis[1] - y_axis[0] if len(y_axis) > 1 else 1.0

    def sample(self, size=1):
        """Sample grid indices exactly, optionally jitter within cells."""
        indices = np.random.choice(len(self.pmf), size=size, p=self.pmf)
        i_x = indices // self.ny
        i_y = indices % self.ny
        
        # Get cell centers
        x_samp = self.x_axis[i_x]
        y_samp = self.y_axis[i_y]
        
        # Optional jittering: uniform within ±0.5*cell
        if self.jitter:
            x_samp = x_samp + np.random.uniform(-0.5*self.dx, 0.5*self.dx, size=size)
            y_samp = y_samp + np.random.uniform(-0.5*self.dy, 0.5*self.dy, size=size)
        
        return np.column_stack([x_samp, y_samp])
