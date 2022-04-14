# Interoperability Between Julia and Python

There are currently two ways to call Julia from Python and vice versa.

## JuliaCall

Vlasiator.jl can be called from Python via [JuliaCall](https://cjdoris.github.io/PythonCall.jl/dev/juliacall/).
JuliaCall will link to the first Julia version in the system path. If Vlasiator.jl has been installed, we can use it directly; otherwise we need to state it in the [juliacalldeps.json](https://cjdoris.github.io/PythonCall.jl/dev/juliacall/#juliacalldeps.json) file.

```python
from juliacall import Main as jl
jl.seval("using Vlasiator")
file = "bulk.1d.vlsv"
meta = jl.load(file)
```

Matplotlib can then be used to visualize the data.

```python
from matplotlib import pyplot as plt
import numpy as np
rho = jl.readvariable(meta, "proton/vg_rho")
x = np.arange(meta.coordmin[0], meta.coordmax[0], meta.dcoord[0])
plt.plot(x, rho)
plt.show()
```

!!! warn
    There seems to be issue that JuliaCall may decide to check pkg installation every time for a new session. We need to first make sure that PythonCall is installed in Julia; then make sure your `PYTHONPATH` is properly set. See this [issue](https://github.com/cjdoris/PythonCall.jl/issues/144) for more information.

## PyJulia

Vlasiator.jl can also be called from Python with the aid of [PyJulia](https://pyjulia.readthedocs.io/en/latest/).
Following the installation steps described in the manual[^1], and then inside Python REPL:

```python
# Handling initialization issue for Conda
from julia.api import Julia
jl = Julia(compiled_modules=False)

from julia import Vlasiator
file = "bulk1.0001000.vlsv"
meta = Vlasiator.load(file)
var = "proton/vg_rho"
data = Vlasiator.readvariable(meta, var)
```

To run a Julia script in Python,

```python
# Handling initialization issue for Conda
from julia.api import Julia
jl = Julia(compiled_modules=False)
jl.eval('include("examples/demo_2dplot_pyplot.jl")')
import matplotlib.pyplot as plt
plt.show()
```

!!! note
    This approach is for you to have a taste of the package with a Python frontend. The workaround shown above for handling the static python libraries makes it slow for regular use. An alternative solution would be creating system images, but as of Julia 1.6 the user experience is not smooth. For better integrated experience with its full power, it is recommended to use the package inside Julia.

[^1]: For Debian-based Linux distributions, it gets a little bit tricky. Please refer to [Troubleshooting](https://pyjulia.readthedocs.io/en/latest/troubleshooting.html) for details.