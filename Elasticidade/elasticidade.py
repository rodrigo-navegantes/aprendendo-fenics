from fenics import *
import matplotlib.pyplot as plt

E = 10.0                                                                            #Módulo de Young
nu = 0.3                                                                            #Coeficiente de Poisson

mu = E / (2.0 * (1.0 + nu))                                                         #Módulo de Cisalhamento
lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))                                    #Parametro de Lamé

mesh = RectangleMesh(Point(0, 0), Point(2, 1), 40, 20)

V = VectorFunctionSpace(mesh, "P", 2)                                               #Desta vez tratamos de um espaço vetorial

# Condições de contorno
left = CompiledSubDomain("near(x[0], 0.0)")                                         #Identifica a Borda Esquerda (Fixa)
right = CompiledSubDomain("near(x[0], 2.0)")                                        #Identifica a Borda Direita (Submetida a tração)                

bc_left = DirichletBC(V, Constant((0.0, 0.0)), left)

traction = Constant((0.1, 0.0))                                             

u = TrialFunction(V)
v = TestFunction(V)

def epsilon(u): return sym(grad(u))                                                 # Tensor de Deformação
def sigma(u): return lmbda * tr(epsilon(u)) * Identity(2) + 2.0 * mu * epsilon(u)   # Tensor de Tensão

# Forma fraca
a = inner(sigma(u), epsilon(v)) * dx
L = dot(traction, v) * ds

# Marcar as fronteiras (?)
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.set_all(0)
right.mark(boundaries, 1)
ds = Measure("ds", domain=mesh, subdomain_data=boundaries)

# Aplicamos o solver
u_sol = Function(V)
solve(a == L, u_sol, bc_left)

# Plotar campo de deslocamento
plt.figure()
plot(u_sol, title="Deslocamento u(x, y)")

u_magnitude = sqrt(dot(u_sol, u_sol))                                               # ||u||
V0 = FunctionSpace(mesh, 'P', 1)                                                    # Definimos outro espaço escalar
u_mag_proj = project(u_magnitude, V0)                                               # Projetamos a norma do vetor no espaço

plt.figure()
plot(u_mag_proj, title="Magnitude do Deslocamento |u|")
plt.show()  
