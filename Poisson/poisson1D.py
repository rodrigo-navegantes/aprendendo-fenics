from fenics import *

# FENICs built-in mesh
mesh = UnitIntervalMesh(32)

# Espaço de funções
V = FunctionSpace(mesh, "Lagrange", 1)  # Funções de Lagrange

# Definir condição de contorno (u = 0 na fronteira)
u_0 = Constant(0.0)
u_1 = Constant(1.0)

left = CompiledSubDomain("near(x[0], 0.0)") #Identifica a Borda Esquerda
right = CompiledSubDomain("near(x[0], 1.0)") #Identifica a Borda Direita

bc_left = DirichletBC(V, u_0, left)
bc_right = DirichletBC(V, u_1, right)

# Definir o problema variacional
u = TrialFunction(V)
v = TestFunction(V)

# Termo fonte
f = Expression(
    "exp(3*x[0])",
    degree=5
) 

# Problema variacional (forma fraca do PDE),
# deve ser escrito na forma: a(u, v) = L(v)
a = dot(grad(u), grad(v)) * dx
L = f * v * dx

# Resolve o problema
u_sol = Function(V)
solve(a == L, u_sol, [bc_left, bc_right])

# Salva a solução em formato .pvd (para abrir no ParaView)
vtkfile = File("output/solution1D.pvd")
vtkfile << u_sol

# Calcula o erro L2
u_exact = Expression(
    "(1 - exp(3*x[0]) +8* x[0] + exp(3)*x[0])/9", 
    degree=5
)
error_L2 = errornorm(u_exact, u_sol, 'L2')
print(f"L2 error: {error_L2:.2e}")


