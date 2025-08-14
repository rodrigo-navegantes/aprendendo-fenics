from dolfin import *
import numpy as np

# Criação da malha e espaços de função
malha = UnitSquareMesh(40, 40)                
V = FunctionSpace(malha, "CG", 1)             
Vmalha = VectorFunctionSpace(malha, "CG", 1)

# Condição de contorno: u = 0 em todo o contorno
u_contorno = Constant(0.0)
condicao_contorno = DirichletBC(V, u_contorno, "on_boundary") # aqui usamos a função pre-definida "on_boundary"

# Problema variacional para Laplace: Δu = -f
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(1.0)
a = dot(grad(u), grad(v)) * dx
L = f * v * dx

u_sol = Function(V)
deslocamento = Function(Vmalha)

arquivo_solucao = File("ALE/output/solucao.pvd")
arquivo_malha = File("ALE/output/malha.pvd")

# Parâmetros do passo de tempo
num_passos = 40
dt = 0.05

# Loop no tempo
for passo in range(num_passos):
    tempo = (passo + 1) * dt

    # Construimos um campo de deslocamento dependente do tempo
    coordenadas = malha.coordinates()     
    vetor_desloc = np.zeros(coordenadas.shape)

    # Exemplo: deslocamento horizontal, maior próximo ao meio em y
    vetor_desloc[:, 0] = 0.005 * np.sin(2.0 * np.pi * coordenadas[:, 1]) * np.sin(2.0 * np.pi * tempo / (num_passos * dt))
    vetor_desloc[:, 1] = 0

    # Atribuir o vetor inteiro de uma vez
    deslocamento.vector()[:] = vetor_desloc.flatten()

    # Movemos a malha usando ALE
    ALE.move(malha, deslocamento)
    malha.bounding_box_tree().build(malha)  # atualizar estrutura de busca

    # Resolver o PDE na malha movida
    solve(a == L, u_sol, condicao_contorno)

    # Salvar solução e malha com valor de tempo
    arquivo_solucao << (u_sol, tempo)
    arquivo_malha << (malha, tempo)

    print(f"Passo {passo+1}/{num_passos} concluído (t = {tempo:.3f})")