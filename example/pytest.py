import numpy as np
from pathlib import Path
from pyszx import SZx
import sys

# prepare your data in numpy array format
HOME = str(Path.home())
data = np.random.rand(10000).astype(np.float32)

# init SZx
# Please change the path to the SZx dynamic library file in your system
lib_extention = {
    "darwin": "libSZx.dylib",
    "windows": "SZx.dll",
}.get(sys.platform, "libSZ3c.so")

szx = SZx("../build/lib/{}".format(lib_extention))

# compress, both input and output data are numpy array
data_cmpr, cmpr_ratio = szx.compress(data, 1, 0, 1e-2)
print("compression ratio = {:5G}".format(cmpr_ratio))

# decompress, both input and output data are numpy array
data_dec = szx.decompress(data_cmpr, data.shape, data.dtype)

# verify
szx.verify(data, data_dec)
