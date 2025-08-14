from fenics import *

t_end = 100.0   #Tempo final de simulação
dt = 0.1        #Time step
k = 3000        #Condutividade Térmica
u_in = 20       #Temperatura na Extremidade Superior
u_out = -20     #Temperatura na Extremidade Inferior

xml_file = "/home/rodrigo/Documents/IC/FEniCS/Helix/helix.xml"
mesh = Mesh(xml_file)
fd = MeshFunction('size_t', mesh, "/home/rodrigo/Documents/IC/FEniCS/Helix/helix_facet_region.xml");

V = FunctionSpace(mesh, 'P', 1)

bc1 = DirichletBC(V, Constant(u_in), fd, 3) 
bc2 = DirichletBC(V, Constant(u_out), fd, 2)
bc = [bc1, bc2]

u = TrialFunction(V)
v = TestFunction(V)
u_n = Function(V)

F = u*v*dx + dt*k*dot(grad(u), grad(v))*dx - u_n*v*dx
a, L = lhs(F), rhs(F)

u = Function(V)
t = 0
vtkfile = File('/home/rodrigo/Documents/IC/FEniCS/Helix/output/output.pvd')

num_steps = int(t_end/dt)
for n in range(num_steps):
    t += dt
    solve(a == L, u, bc)
    u_n.assign(u)
    vtkfile << (u, t)
