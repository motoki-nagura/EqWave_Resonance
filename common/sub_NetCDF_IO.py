import numpy as np
from netCDF4 import Dataset

def infer_dim_names(shape, coord_info=None):
    """
    Infer dimension names from shape with optional user hints.
    """
    default_order = ["time", "level", "lat", "lon"]

    ndim = len(shape)
    dim_names = []

    for i, size in enumerate(shape):
        if coord_info and i in coord_info:
            dim_names.append(coord_info[i]["name"])
        else:
            # Heuristic mapping (last dims usually spatial)
            if ndim <= len(default_order):
                dim_names.append(default_order[-ndim + i])
            else:
                dim_names.append(f"dim_{i}")

    return dim_names


def create_or_get_dim(ncfile, name, size):
    if name not in ncfile.dimensions:
        ncfile.createDimension(name, size)


def create_coordinate_variable(ncfile, name, size, values=None):
    if name in ncfile.variables:
        return ncfile.variables[name]

    var = ncfile.createVariable(name, "f8", (name,))

    if values is not None:
        var[:] = values
    else:
        var[:] = np.arange(size)

    # Minimal CF-style metadata
    if name == "time":
        var.units = "index"
    elif name == "lat":
        var.units = "degrees_north"
    elif name == "lon":
        var.units = "degrees_east"

    return var


def write_variable_to_netcdf(
    array,
    filename,
    var_name,
    coord_info=None,
    mode="a"
):
    """
    Write or append a variable with automatic dimension + coordinate handling.

    Parameters:
        array (np.ndarray)
        filename (str)
        var_name (str)
        coord_info (dict):
            Optional per-dimension metadata, e.g.:
            {
                0: {"name": "time", "values": np.array([...])},
                1: {"name": "lat", "values": lat_array},
                2: {"name": "lon", "values": lon_array}
            }
        mode (str): "w" (new file) or "a" (append)
    """

    if not isinstance(array, np.ndarray):
        raise TypeError("Input must be a NumPy ndarray")

    with Dataset(filename, mode, format="NETCDF4") as ncfile:

        shape = array.shape
        dim_names = infer_dim_names(shape, coord_info)

        # Create dimensions + coordinate variables
        for i, (dim_name, size) in enumerate(zip(dim_names, shape)):
            create_or_get_dim(ncfile, dim_name, size)

            coord_values = None
            if coord_info and i in coord_info:
                coord_values = coord_info[i].get("values")

            create_coordinate_variable(ncfile, dim_name, size, coord_values)

        # Avoid overwriting existing variable
        if var_name in ncfile.variables:
            raise ValueError(f"Variable '{var_name}' already exists")

        # Create variable
        var = ncfile.createVariable(
            var_name,
            array.dtype,
            tuple(dim_names),
            zlib=True,
            complevel=4
        )

        var[:] = array

        # Attach coordinate reference
        var.coordinates = " ".join(dim_names)

        print(f"OUT >> {filename}:  Variable '{var_name}' with dims {dim_names}")


def read_variable_from_netcdf(
    filename,
    var_name,
    squeeze=True,
    return_coords=False,
    slices=None
):
    """
    Read a specific variable from a NetCDF file.

    Parameters:
        filename (str): Path to NetCDF file
        var_name (str): Variable name to read
        squeeze (bool): Remove singleton dimensions
        return_coords (bool): Also return coordinate variables
        slices (tuple or None): Optional slicing (e.g., (0, :, :) )

    Returns:
        data (np.ndarray)
        coords (dict, optional)
    """

    with Dataset(filename, "r") as ncfile:

        if var_name not in ncfile.variables:
            raise KeyError(f"Variable '{var_name}' not found in file")

        var = ncfile.variables[var_name]

        # Read data (with optional slicing)
        if slices is not None:
            data = var[slices]
        else:
            data = var[:]

        # Convert masked array → ndarray (optional)
        if hasattr(data, "mask"):
            data = np.array(data.filled(np.nan))

        if squeeze:
            data = np.squeeze(data)

        if not return_coords:
            print(f"IN << {filename}:  Variable '{var_name}'")
            return data

        # Extract coordinate variables
        coords = {}
        for dim in var.dimensions:
            if dim in ncfile.variables:
                coords[dim] = ncfile.variables[dim][:]

        print(f"IN << {filename}:  Variable '{var_name}' with dims {var.dimensions}")

        return data, coords


def main_test():
    """
       Example Usage
    """
    # 1. Create File + First Variable

    temp = np.random.rand(10, 50, 100)  # time, lat, lon

    lat = np.linspace(-90, 90, 50)
    lon = np.linspace(0, 360, 100)

    write_variable_to_netcdf(
        temp,
        "climate.nc",
        "temperature",
        coord_info={
            0: {"name": "time"},
            1: {"name": "lat", "values": lat},
            2: {"name": "lon", "values": lon},
        },
        mode="w"
    )

    # 2. Append Another Variable

    pressure = np.random.rand(10, 50, 100)

    write_variable_to_netcdf(
        pressure,
        "climate.nc",
        "pressure",
        coord_info={
            0: {"name": "time"},
            1: {"name": "lat"},
            2: {"name": "lon"},
        },
        mode="a"
    )

    # 3. Simple Read

    data = read_variable_from_netcdf("climate.nc", "temperature")
    print('  data.shape = ',data.shape)

    # 4. Read With Coordinates

    data, coords = read_variable_from_netcdf(
        "climate.nc",
        "temperature",
        return_coords=True
    )

    lat = coords["lat"]
    lon = coords["lon"]
    time = coords["time"]

    print('  coords.keys() = ',coords.keys())  # e.g., dict_keys(['time', 'lat', 'lon'])

    # 5. Partial Read (Efficient for Large Files)

    # Read only first timestep
    data = read_variable_from_netcdf(
        "climate.nc",
        "temperature",
        slices=(0, slice(None), slice(None))
    )
    print('  data.shape = ',data.shape)


if __name__ == "__main__":
    import sys
    sys.exit(main_test())

