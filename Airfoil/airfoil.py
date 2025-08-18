from dolfin import *

# -----------------------------
# CARREGAR MALHA DO AEROFÓLIO
# -----------------------------
mesh = Mesh()
with XDMFFile("/home/rodrigo/Documents/IC/FEniCS/Airfoil/malha/airfoil.xdmf") as infile:
    infile.read(mesh)

# Carregar marcadores de fronteira (Physical Curves)
mvc = MeshValueCollection("size_t", mesh, 1)
with XDMFFile("/home/rodrigo/Documents/IC/FEniCS/Airfoil/malha/airfoil.xdmf") as infile:
    infile.read(mvc, "Aerofolio")  # "Aerofolio" é o nome do Physical Curve
boundaries = MeshFunction("size_t", mesh, mvc)

# -----------------------------
# ESPAÇOS DE FUNÇÃO (P2-P1 Taylor-Hood)
# -----------------------------
# Criar elementos separados e depois combiná-los
V_el = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # Velocidade (P2)
Q_el = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # Pressão (P1)
W = FunctionSpace(mesh, MixedElement([V_el, Q_el]))   # Espaço misto

# Funções de teste e solução
(v, q) = TestFunctions(W)
w = Function(W)
(u, p) = split(w)

# -----------------------------
# CONDIÇÕES DE CONTORNO
# -----------------------------
bcs = []

# No-slip no aerofólio (ID = 1)
bcs.append(DirichletBC(W.sub(0), Constant((0,0)), boundaries, 1))

# Velocidade de entrada (ID = 103)
bcs.append(DirichletBC(W.sub(0), Constant((1.0,0.0)), boundaries, 103))

# Topo e base podem ser livres ou como paredes, opcional
# bcs.append(DirichletBC(W.sub(0), Constant((0,0)), boundaries, 102))  # topo
# bcs.append(DirichletBC(W.sub(0), Constant((0,0)), boundaries, 100))  # base

# -----------------------------
# PARÂMETROS FÍSICOS
# -----------------------------
mu = 1.0
rho = 1.0
g = Constant((0, -9.81))

# -----------------------------
# FORMULAÇÃO FRACA DE STOKES
# -----------------------------
epsilon = lambda u: sym(grad(u))
sigma   = lambda u, p: 2*mu*epsilon(u) - p*Identity(len(u))

F = (inner(sigma(u, p), grad(v)) - rho*dot(g, v) - div(u)*q)*dx

# -----------------------------
# RESOLVER O SISTEMA
# -----------------------------
solve(F == 0, w, bcs)

# Separar velocidade e pressão
(u, p) = w.split()

# -----------------------------
# SALVAR PARA PARAVIEW
# -----------------------------
File("/home/rodrigo/Documents/IC/FEniCS/Airfoil/output/velocidade.pvd") << u
File("/home/rodrigo/Documents/IC/FEniCS/Airfoil/output/pressao.pvd") << p