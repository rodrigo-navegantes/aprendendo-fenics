import numpy as np

# Load coordinates, skipping the first line with "NACA 4412"
data = np.loadtxt("airfoil.dat", skiprows=1)

# Remove duplicate trailing edge point at the end
if np.allclose(data[0], data[-1]):
    data = data[:-1]

# Write Gmsh .geo file
with open("airfoil.geo", "w") as f:
    # Points
    for i, (x, y) in enumerate(data, start=1):
        f.write(f"Point({i}) = {{{x:.6f}, {y:.6f}, 0}};\n")
    
    # Closed line
    n_points = len(data)
    line_str = "Line(1) = {" + ", ".join(str(i) for i in range(1, n_points + 1)) + ", 1};\n"
    f.write(line_str)

