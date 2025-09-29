from fenics import *

# Criação da malha e espaço de funções
mesh = UnitSquareMesh(32, 32) #esta opção cria sempre triângulos
V = FunctionSpace(mesh, "Lagrange", 1)  # Funções de Lagrange

# Definir condição de contorno (u = 0 na fronteira)
u_D = Constant(0.0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Definir o problema variacional
u = TrialFunction(V)
v = TestFunction(V)
f = Expression(
    "10*exp(-5*((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5)))",
    degree=5
)  # Função gaussiana centrada em (0,5 0,5) com pico em 10

L = f * v * dx
a = dot(grad(u), grad(v)) * dx

# Resolver o problema
u_sol = Function(V)
solve(a == L, u_sol, bc)

# Salvar em formato .pvd (para abrir no ParaView)
vtkfile = File("output/solution.pvd")
vtkfile << u_sol
