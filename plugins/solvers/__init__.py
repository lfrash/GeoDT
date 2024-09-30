# Define the __all__ variable
__all__ = ["geodt", "gringarten"]

# Import the submodules
# from .libs.gt_linalg import solve
# from .libs import gt_vtk as sg
# from .libs import gt_properties
# water = gt_properties.water()
from . import geodt
from . import gringarten
from . import libs
