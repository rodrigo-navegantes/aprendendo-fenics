from fenics import *
import os

# Disable HDF5 locking (helps on some systems)
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

# --- Load mesh ---
mesh = Mesh()
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)

# --- Load boundaries ---
mvc = MeshValueCollection("size_t", mesh, 1)  # 1 = facets
with XDMFFile("mf.xdmf") as infile:
    infile.read(mvc, "gmsh:physical")
boundaries = MeshFunction("size_t", mesh, mvc)

# --- Function space ---
V = FunctionSpace(mesh, "Lagrange", 1)

# --- Dirichlet BCs (IDs match .geo physical tags) ---
bc_left  = DirichletBC(V, Constant(0.0), boundaries, 4)  # left boundary
bc_right = DirichletBC(V, Constant(0.0), boundaries, 2)  # right boundary
bc_bottom = DirichletBC(V, Constant(0.0), boundaries, 5)  # bottom boundary
bc_top = DirichletBC(V, Constant(0.0), boundaries, 3)  # top boundary

# --- Variational problem ---
u = TrialFunction(V)
v = TestFunction(V)
f = Expression(
    "(x[0]-1)", degree=5
)

a = dot(grad(u), grad(v)) * dx
L = f * v * dx

# --- Solve ---
u_sol = Function(V)
solve(a == L, u_sol, [bc_left, bc_right, bc_bottom, bc_top])

# --- Save solution ---
os.makedirs("output", exist_ok=True)
File("output/solution.pvd") << u_sol
