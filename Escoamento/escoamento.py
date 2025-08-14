from fenics import *

nu = 1.0  # viscosidade

# Malha
mesh = RectangleMesh(Point(0, 0), Point(2, 1), 64, 32)

V = VectorElement("P", mesh.ufl_cell(), 2)   # Espaço vetorial para velocidade
Q = FiniteElement("P", mesh.ufl_cell(), 1)   # Espaço escalar para pressão
W = FunctionSpace(mesh, MixedElement([V, Q]))  # Espaço misto W = V * Q

# Funções de teste e tentativa
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

# Força externa
f = Constant((0.0, 0.0))

# Condições de contorno
inflow = Expression(("4.0 * x[1] * (1.0 - x[1])", "0.0"), degree=2)  # perfil parabólico
noslip = Constant((0.0, 0.0))                                       # não deslizamento

def left(x, on_boundary):  return near(x[0], 0) and on_boundary
def right(x, on_boundary): return near(x[0], 2) and on_boundary
def walls(x, on_boundary): return (near(x[1], 0) or near(x[1], 1)) and on_boundary

bc_in   = DirichletBC(W.sub(0), inflow, left)
bc_wall = DirichletBC(W.sub(0), noslip, walls)
bc_out  = DirichletBC(W.sub(1), Constant(0), right)  # pressão = 0 na saída

bcs = [bc_in, bc_wall, bc_out]

# Formulação variacional (Stokes)
a = (nu * inner(grad(u), grad(v)) - div(v) * p + q * div(u)) * dx
L = inner(f, v) * dx

# Resolver o sistema
w = Function(W)
solve(a == L, w, bcs)

# Separar velocidade e pressão
u_sol, p_sol = w.split()

# Salvar resultados para o ParaView
File("Escoamento/output/velocity.pvd") << u_sol
File("Escoamento/output/pressure.pvd") << p_sol
