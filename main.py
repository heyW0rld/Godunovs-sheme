#ToDo
#Создавать слой, задавая количество ячеек в нем в конструкторе

import matplotlib.pyplot as plt
import math
import copy
from decimal import Decimal
	


from layer import Layer
from cell import Cell
from layer import Layer
from pur import PUR
import calc_pur

layers = []                 #Расчетная сетка

l = Decimal(0.1)            #начальное положение
N = 10                      #количество разбиений
h_0 = l / N                 #шаг по кооординате на 0 слое
tau_0 = Decimal(10e-4)      #шаг по времени на 0 слое
kappa = Decimal(1.4)        #показатель адиабаты
S = Decimal(0.03)           #сечение снаряда
M = Decimal(3)              #масса снаряда
L = Decimal(3)              #длина ствола


#Задание начальных параметров
def init():
    layer = Layer()

    layer.h = h_0
    layer.tau = tau_0
    layer.x = l

    current_x = Decimal(0)
    for i in range(N + 1):    
        cell = Cell()
        pur_midl = PUR()
        pur = PUR()
        cell.x = current_x
        current_x += h_0
        pur.p, pur.u, pur.r = Decimal(1), Decimal(0), Decimal(1)
        pur_midl.p, pur_midl.u, pur_midl.r = Decimal(1), Decimal(0), Decimal(1)
        cell.pur_midl = pur
        cell.pur = pur
        layer.append_cell(cell)

    layers.append(layer)
       
        
#Вычисление больших величин
def calc_PUR(cell_prev: Cell, cell_current: Cell, cell_prev_next: Cell, cell_current_next: Cell):

    #Вычисление координаты на следующем слое
    # def x_n1() -> Decimal:
    #     def c(p, r) -> Decimal:
    #         return Decimal.sqrt(kappa * p / r)

    #     def p() -> Decimal:
    #         p_prev = cell_prev.pur_midl.p
    #         r_prev = cell_prev.pur_midl.r
    #         u_cur = cell_current.pur.u
    #         u_prev = cell_prev.pur_midl.u
    #         return  p_prev * (1 - (kappa - 1) / (2 * c(p_prev, r_prev)) * (u_cur - u_prev)) ** (2 * kappa / (kappa - 1))

    #     def ro() -> Decimal:
    #         p_prev = cell_prev.pur_midl.p
    #         r_prev = cell_prev.pur_midl.r
    #         return r_prev * (p() / p_prev) ** (1 / kappa)

    #     x_cur = cell_current.x
    #     u_cur = cell_current.pur.u   
    #     tau = layers[-2].tau

    #     return x_cur + (u_cur + 2 * c(p(), ro()) / (kappa - 1)) * tau + 2 * c(p(), ro()) ** 2 / (kappa-1) * M / (S * p())* \
    #                                                         (1 - (1 + (kappa + 1) / (2 * c(p(), ro())) * p() * S / M * tau) ** (2 / (kappa + 1)))
    
    u_m, p_m = calc_pur.formul_1_2(cell_prev, cell_current)
    ro_m = Decimal(0)
    cell_current.pur.p = p_m
    cell_current.pur.u = u_m

    # cell_current_next.x = x_n1()
    
    deltax = cell_current_next.x - cell_current.x

    p_prev = cell_prev.pur_midl.p
    p_next = cell_current.pur_midl.p
    r_prev = cell_prev.pur_midl.r
    r_next = cell_current.pur_midl.r
    tau = layers[-2].tau
    u_cur = cell_current.pur.u

    if p_m <= p_prev and p_m <= p_next:
        if deltax/tau > u_cur:
            ro_m = r_prev*(p_m/p_prev)**(1/kappa)
        else:
            ro_m = r_next*(p_m/p_next)**(1/kappa)

    elif p_m > p_prev and p_m <= p_next:
        u_m, p_m = calc_pur.formul_1_4(cell_prev, cell_current)
        if deltax/tau > u_cur:
            ro_m = r_prev*((kappa+1)*p_m + (kappa-1)*p_prev) / \
                ((kappa - 1)*p_m + (kappa+1)*p_prev)
        else:
            ro_m = r_next*(p_m/p_next)**(1/kappa)

    elif p_m <= p_prev and p_m > p_next:
        u_m, p_m = calc_pur.formul_1_3(cell_prev, cell_current)
        if deltax/tau > u_cur:
            ro_m = r_prev*(p_m/p_prev)**(1/kappa)
        else:
            ro_m = r_next*((kappa+1)*p_m + (kappa-1)*p_next) / \
                ((kappa - 1)*p_m + (kappa+1)*p_next)

    elif p_m > p_prev and p_m > p_next:
        u_m, p_m = calc_pur.formul_3_4(cell_prev, cell_current)
        if deltax/tau > u_cur:
            ro_m = r_prev*((kappa+1)*p_m + (kappa-1)*p_prev) / \
                ((kappa - 1)*p_m + (kappa+1)*p_prev)
        else:
            ro_m = r_next*((kappa+1)*p_m + (kappa-1)*p_next) / \
                ((kappa - 1)*p_m + (kappa+1)*p_next)

    cell_current.pur.p = p_m
    cell_current.pur.u = u_m
    cell_current.pur.r = ro_m


#Вычисление последней ячейки
def calc_last_cell(cell_prev: Cell, cell_current: Cell, cell_prev_next: Cell, cell_current_next: Cell):
    x_cur = cell_current.x
    u_prev = cell_prev.pur_midl.u
    u_cur = cell_current.pur.u 
    p_prev = cell_prev.pur_midl.p
    r_prev = cell_prev.pur_midl.r  
    tau = layers[-1].tau

    def c(p, r) -> Decimal:
        return Decimal.sqrt(kappa * p / r)

    def p() -> Decimal:
        return  p_prev * (1 - (kappa - 1) / (2 * c(p_prev, r_prev)) * (u_cur - u_prev)) ** (2 * kappa / (kappa - 1))

    def ro() -> Decimal:
        return r_prev * (p() / p_prev) ** (1 / kappa)

    def x_n1() -> Decimal:
        return x_cur + (u_cur + 2 * c(p(), ro()) / (kappa - 1)) * tau + 2 * c(p(), ro()) ** 2 / (kappa-1) * M / (S * p())* \
                                                            (1 - (1 + (kappa + 1) / (2 * c(p(), ro())) * p() * S / M * tau) ** (2 / (kappa + 1)))

    def u_n1() -> Decimal:
        return u_cur + 2 * c(p(), ro()) / (kappa - 1) * (1 - (1 + (kappa + 1) / (2 * c(p(), ro())) * p() * S / M * tau) \
                                                                                        ** (-(kappa - 1) / (kappa + 1)))
    
    cell_current_next.x = x_n1()
    
    pur = PUR()
    pur.u = u_n1()

    p_prev_n1 = cell_prev_next.pur_midl.p
    u_cur_n1 = pur.u
    u_prev_n1 = cell_prev_next.pur_midl.u   
    def p_n1() -> Decimal:
        return p_prev_n1 * (1 - (kappa - 1) / (2 * c(p_prev, r_prev)) * (u_cur_n1 - u_prev_n1)) ** (2 * kappa / (kappa - 1))
    pur.p = p_n1()

    r_prev_n1 = cell_prev_next.pur_midl.r
    p_cur_n1 = pur.p
    def r_n1() -> Decimal:
        return  r_prev_n1 * (p_cur_n1 / p_prev_n1) ** (1/kappa)

    pur.r = r_n1()      

    cell_current_next.pur = pur


#Вычисление промежуточных значений на следующем слое
def calc_next_midl_cell():
    cur_layer = layers[-2]
    next_layer = layers[-1]
    tau_cur = layers[-2].tau

    for i in range(1, N):
        r_midl = cur_layer.cells[i].pur_midl.r
        r_next = cur_layer.cells[i+1].pur.r
        r_cur = cur_layer.cells[i].pur.r
        h = cur_layer.h
        h1 = next_layer.h
        u_next = cur_layer.cells[i+1].pur.u
        u_cur = cur_layer.cells[i].pur.u
        p_midl = cur_layer.cells[i].pur_midl.p
        p_next = cur_layer.cells[i+1].pur.p
        p_cur = cur_layer.cells[i].pur.p
        delta_x_next = next_layer.cells[i+1].x - cur_layer.cells[i+1].x
        delta_x_cur = next_layer.cells[i].x - cur_layer.cells[i].x
        ru_next = cur_layer.cells[i+1].pur.r * cur_layer.cells[i+1].pur.u
        ru_cur = cur_layer.cells[i].pur.r * cur_layer.cells[i].pur.u

        def r_next_midl() -> Decimal:
            return r_midl * h / h1 + r_next * delta_x_next / h1 - r_cur * delta_x_cur / h1 - (ru_next - ru_cur) * tau_cur / h1

        u_midl = cur_layer.cells[i].pur_midl.u

        def ru_next_midl() -> Decimal:
            return r_midl * u_midl * h / h1 + ru_next * delta_x_next / h1 - ru_cur * delta_x_cur / h1 - ((r_next * u_next ** 2 + p_next) -  \
                                                                                                            (r_cur * u_cur ** 2 + p_cur)) * tau_cur / h1 

        def e(p, r) -> Decimal:
            return p / (r * (kappa - 1))

        def reu2() -> Decimal:
            return r_midl * (e(p_midl, r_midl) + u_midl ** 2 / 2) * h / h1 + (r_next * (e(p_next, r_next) + u_next ** 2 / 2)) * delta_x_next / h1 - \
                                                                                        (r_cur * (e(p_cur, r_cur) + u_cur ** 2 / 2)) * delta_x_cur / h1 - \
                                                                                            ((r_next * u_next * (e(p_next, r_next) + u_next ** 2 / 2) + p_next * u_next) - \
                                                                                                (r_cur * u_cur * (e(p_cur, r_cur) + u_cur ** 2 / 2) + p_cur * u_cur)) * tau_cur / h1

        def u_next_midl() -> Decimal:
            return ru_next_midl() / r_next_midl()

        def e_next_midl() -> Decimal:
            return reu2() / r_next_midl() - u_next_midl() ** 2 / 2

        def p_next_midl() -> Decimal:
            return (kappa - 1) * r_next_midl() * e_next_midl()

        
        pur_midl = PUR()
        pur_midl.p = p_next_midl()
        pur_midl.u = u_next_midl()
        pur_midl.r = r_next_midl()
        next_layer.cells[i].pur_midl = pur_midl 
        
def x1(cell_prev, cell_current) -> Decimal:
    def c(p, r) -> Decimal:
        return Decimal.sqrt(kappa * p / r)
    def p() -> Decimal:
        return  p_prev * (1 - (kappa - 1) / (2 * c(p_prev, r_prev)) * (u_cur - u_prev)) ** (2 * kappa / (kappa - 1))
    def ro() -> Decimal:
        return r_prev * (p() / p_prev) ** (1 / kappa)
        


    x_cur = cell_current.x
    u_prev = cell_prev.pur_midl.u
    u_cur = cell_current.pur.u 
    p_prev = cell_prev.pur_midl.p
    r_prev = cell_prev.pur_midl.r  
    tau = layers[-1].tau

    return x_cur + (u_cur + 2 * c(p(), ro()) / (kappa - 1)) * tau + 2 * c(p(), ro()) ** 2 / (kappa-1) * M / (S * p())* \
                                                            (1 - (1 + (kappa + 1) / (2 * c(p(), ro())) * p() * S / M * tau) ** (2 / (kappa + 1)))

#Основной алгоритм

y = [0.1]
x = [0]
p = [1]
u = [0]
r = [1]

def main():
    #Создание и заполнение первого слоя
    init()
    n = 0

    #Цикл заполнения ячеек в слое 
    while layers[-1].x < L:
        #Создание нового слоя
        new_layer = Layer()
        new_layer.tau = tau_0
        for i in range(N+1):
            new_layer.append_cell(Cell())
        layers.append(new_layer)

        #Значение координаты на правой границе
        layers[-1].cells[N].x = x1(layers[-2].cells[N - 1], layers[-2].cells[N])
        layers[-1].x = layers[-1].cells[N].x
        new_layer.h = layers[-1].cells[N].x / N
        current_x = 0
        for i in range(N+1):
            layers[-1].cells[i].x = current_x
            current_x += layers[-1].h

        #Находим большие величины P U R
        #Установка половинных PUR для первой ячейки
        nullPUR = PUR()
        firstCell = layers[-2].cells[1] 
        nullPUR.p = firstCell.pur_midl.p
        nullPUR.u = -firstCell.pur_midl.u
        nullPUR.r = firstCell.pur_midl.r
        layers[-2].cells[0].pur_midl = nullPUR

        #От 0 до N
        for i in range(1, N):
            calc_PUR(layers[-2].cells[i-1], layers[-2].cells[i], layers[-1].cells[i-1], layers[-1].cells[i]) 
        
        #N-я ячейка 
        if(n == 0):
            prelast_PUR = copy.copy(layers[-2].cells[N-1].pur_midl)
            layers[-2].cells[N].pur = prelast_PUR 


        calc_next_midl_cell()
        calc_last_cell(layers[-2].cells[N - 1], layers[-2].cells[N], layers[-1].cells[N - 1], layers[-1].cells[N])
        
        n+=1
        print(layers[-1].x)
        y.append(layers[-1].x)
        x.append(y[-1]+tau_0)
        p.append(layers[-1].cells[-1].pur.p)
        u.append(layers[-1].cells[-1].pur.u)
        r.append(layers[-1].cells[-1].pur.r)

     

main()

plt.plot(x,y)
plt.plot(x,p)
plt.plot(x,u)
plt.plot(x,r)
plt.show()