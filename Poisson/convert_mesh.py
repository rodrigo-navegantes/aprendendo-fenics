import meshio

# Read gmsh mesh
msh = meshio.read("square.msh")

# --- Bulk mesh (triangles only) ---
triangle_mesh = meshio.Mesh(
    points=msh.points,
    cells=[("triangle", msh.cells_dict["triangle"])],
    cell_data={"gmsh:physical": [msh.cell_data_dict["gmsh:physical"]["triangle"]]},
)
triangle_mesh.write("mesh.xdmf")

# --- Boundary mesh (lines only) ---
line_mesh = meshio.Mesh(
    points=msh.points,
    cells=[("line", msh.cells_dict["line"])],
    cell_data={"gmsh:physical": [msh.cell_data_dict["gmsh:physical"]["line"]]},
)
line_mesh.write("mf.xdmf")
