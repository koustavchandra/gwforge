from scipy.interpolate import RegularGridInterpolator
import numpy as np
from bilby.core.prior import Prior

# Generic 2D interpolated prior class
class Interped2D(Prior):
    """
    A 2D interpolated Prior that samples from a joint PDF over two variables.
    
    Parameters
    ----------
    name : str
        Name of the prior (e.g., 'mass_pair')
    x_axis : array-like
        1D array of values for the first variable
    y_axis : array-like
        1D array of values for the second variable
    pdf_grid : 2D array
        2D grid of probability densities with shape (len(x_axis), len(y_axis))
    x_name : str, optional
        Name of first variable (for labels, default 'x')
    y_name : str, optional
        Name of second variable (for labels, default 'y')
    """
    
    def __init__(self, name, x_axis, y_axis, pdf_grid, x_name='x', y_name='y'):
        super().__init__(name=name, latex_label=f"${name}$", unit=None)
        self.x_axis = np.asarray(x_axis)
        self.y_axis = np.asarray(y_axis)
        self.x_name = x_name
        self.y_name = y_name
        
        # Normalize
        self.pdf_grid = np.asarray(pdf_grid) / np.sum(pdf_grid)
        
        # Bounds
        self.x_min, self.x_max = float(self.x_axis[0]), float(self.x_axis[-1])
        self.y_min, self.y_max = float(self.y_axis[0]), float(self.y_axis[-1])
        
        # Create interpolator
        self.interpolator = RegularGridInterpolator(
            (self.x_axis, self.y_axis), 
            self.pdf_grid, 
            bounds_error=False, 
            fill_value=0.0
        )
        
        # Max density for rejection sampling
        self.max_pdf = float(np.max(self.pdf_grid))
    
    def prob(self, val):
        """Probability density at (x, y)"""
        val_arr = np.asarray(val).reshape(1, 2)
        return float(self.interpolator(val_arr)[0])
    
    def sample(self, size=1):
        """
        Vectorized rejection sampling

        Returns
        -------
        samples : ndarray
            Array of shape (size, 2) with samples drawn from the 2D distribution
        """
        samples = []
        batch = max(1000, size * 10)
        
        while len(samples) < size:
            # Propose candidates uniformly in bounding box
            x_prop = np.random.uniform(self.x_min, self.x_max, size=batch)
            y_prop = np.random.uniform(self.y_min, self.y_max, size=batch)
            pts = np.column_stack([x_prop, y_prop])
            
            # Evaluate density and do rejection
            p_vals = self.interpolator(pts)
            u = np.random.uniform(0, self.max_pdf, size=batch)
            accept = u < p_vals
            
            if np.any(accept):
                accepted_pts = pts[accept]
                for a in accepted_pts:
                    samples.append(a)
                    if len(samples) >= size:
                        break
        
        return np.asarray(samples)[:size]
