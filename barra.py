"""
Neste exemplo resolveremos a equação de calor: u_t = u_xx + f

em uma barra tal que x ∈ [0.0, 1.0]

com Condições de Contorno de Dirichlet Homogêneas ( u(0,0) = u(1,0) = 0 )

e uma fonte de calor f

com uma Condição Inicial u(x,0) = sen(pi * x)
"""

from fenics import *
import matplotlib.pyplot as plt

#Começamos criando uma malha de 32 elementos no intervalo [0,1]
n_elementos = 32
malha = UnitIntervalMesh(n_elementos)

#Criamos um espaço de funções de elmentos finitos de Lagrange de 1ª ordem
lagrange_space= FunctionSpace(malha, "Lagrange", 1)

#Definindo o valor da função nas extremidades
u_na_extremidade = Constant(0.0)

#A função a seguir retorna se estamos ou não na extremidade
def naFronteira(x, na_extremidade): return na_extremidade

#Implementamos a Condição de Contorno e a Condição Inicial u(x,0) = sen(pi * x)
cc = DirichletBC(lagrange_space, u_na_extremidade, naFronteira)
ci = Expression( "sin(3.141 * x[0])", degree=1)
u_velha = interpolate(ci, lagrange_space) # intepola no espaço de funções (Expression > Function)

#Plotamos a em t = 0.0
plt.figure()
plot(u_velha, label="t=0.0")

#Usamos o método de euler impllícito para discretizarmos no tempo 
step_tempo = 0.1
n_step_tempo = 5

#Problema homogêneo...
fonte = Constant(0.0)

#Definimos as funções de teste (v) e de tentativa (u)
u_tent = TrialFunction(lagrange_space)
v_test = TestFunction(lagrange_space)

forma_fraca_res = (
    u_tent * v_test * dx
    +
    step_tempo * dot(
        grad(u_tent),
        grad(v_test),
    ) * dx
    -
    (
        u_velha * v_test * dx
        +
        step_tempo * fonte * v_test * dx
    )
)

#Separamos nossa EDP em lado esquerdo (lhs) e lado direito (rhs)
forma_fraca_lhs = lhs(forma_fraca_res)
forma_fraca_rhs = rhs(forma_fraca_res)

#u_sol será a função que resolveremos para cada ponto no tempo 
u_sol = Function(lagrange_space)

tempo_atual = 0.0
for i in range(n_step_tempo):
    tempo_atual += step_tempo 

    # Finalmente aplicamos o solver... 
    solve(forma_fraca_lhs == forma_fraca_rhs, u_sol, cc)

    u_velha.assign(u_sol)

    plot(u_sol, label=f"t={tempo_atual:1.1f}")


#Plotamos tudo
plt.legend()
plt.title("Condução de Calor em um Barra \n com Condições de Contorno de Dirichlet Homogêneas")
plt.xlabel("Posição")
plt.ylabel("Temperatura")
plt.grid()
plt.show()
