import math

from layer import Layer
from cell import Cell
from layer import Layer
from pur import PUR
import calc_pur

l = 0.1         #начальное положение
N = 100         #количество разбиений
h_0 = l / N     #шаг по кооординате на 0 слое
tau_0 = 10e-4   #шаг по времени на 0 слое
kappa = 1.4     #показатель адиабаты
S = 0.03        #сечение снаряда
M = 3           #масса снаряда
L = 3           #длина ствола

layers = []     #Расчетная сетка

#Шаг 0 (задание начальных параметров)
def init():
    layer = Layer()

    layer.h = l / N
    layer.tau = tau_0
    layers.x = l

    current_x = 0
    for i in range(N + 1):    
        cell = Cell()
        pur = PUR()

        pur.p = 1
        pur.u = 0
        pur.r = 1
        cell.pur_midl = pur

        cell.x = current_x
        current_x += layer.h

        layer.append_cell(cell)

    layers.append(layer)
       
        
#Шаг 1 (нахождение больших величин)
def calc_PUR(cell_prev: Cell, cell_current: Cell):

    #Вычисление координаты на следующем слое
    def x_n1():
        def c(p, r):
            return math.sqrt(kappa * p / r)

        def p():
            p_prev = cell_prev.pur_midl.p
            r_prev = cell_prev.pur_midl.r
            u_cur = cell_current.pur.u
            u_prev = cell_prev.pur_midl.u
            return  p_prev * (1 - (kappa - 1) / (2 * c(p_prev, r_prev)) * (u_cur - u_prev)) ** (2 * kappa / (kappa - 1))

        def ro():
            p_prev = cell_prev.pur_midl.p
            r_prev = cell_prev.pur_midl.r
            return r_prev * (p() / p_prev) ** (1 / kappa)

        x_cur = cell_current.x
        u_cur = cell_current.pur.u   
        tau = layers[-1].tau

        return x_cur + (u_cur + 2 * c(p(), ro()) / (kappa - 1)) * tau + 2 * c(p(), ro()) ** 2 / (kappa-1) * M / (S * p())* \
                                                            (1 - (1 + (kappa + 1) / (2 * c(p(), ro())) * p() * S / M * tau) ** (2 / (kappa + 1)))
    
    u_m, p_m = calc_pur.formul_1_2(cell_prev, cell_current)
    ro_m = 0
    cell_current.pur.p = p_m
    cell_current.pur.u = u_m

    cell_current.x1 = x_n1()
    
    deltax = cell_current.x1 - cell_current.x

    p_prev = cell_prev.pur_midl.p
    p_next = cell_current.pur_midl.p
    tau = layers[-1].tau
    u_cur = cell_current.pur.u 

    if p_m < p_prev and p_m < p_next:
        if deltax/tau > u_cur:
            ro_m = p_prev*(p_m/p_prev)**(1/kappa)
        else:
            ro_m = p_next*(p_m/p_next)**(1/kappa)

    elif p_m > p_prev and p_m < p_next:
        u_m, p_m = calc_pur.formul_1_4(cell_prev, cell_current)
        if deltax/tau > u_cur:
            ro_m = p_prev*((kappa+1)*p_m + (kappa-1)*p_prev) / \
                ((kappa - 1)*p_m + (kappa+1)*p_prev)
        else:
            ro_m = p_next*(p_m/p_next)**(1/kappa)

    elif p_m < p_prev and p_m > p_next:
        u_m, p_m = calc_pur.formul_1_3(cell_prev, cell_current)
        if deltax/tau > u_cur:
            ro_m = p_prev*(p_m/p_prev)**(1/kappa)
        else:
            ro_m = p_next*((kappa+1)*p_m + (kappa-1)*p_next) / \
                ((kappa - 1)*p_m + (kappa+1)*p_next)

    elif p_m > p_prev and p_m > p_next:
        if deltax/tau > u_cur:
            ro_m = p_prev*((kappa+1)*p_m + (kappa-1)*p_prev) / \
                ((kappa - 1)*p_m + (kappa+1)*p_prev)
        else:
            ro_m = p_next*((kappa+1)*p_m + (kappa-1)*p_next) / \
                ((kappa - 1)*p_m + (kappa+1)*p_next)

    cell_current.pur.r = ro_m


#Вычисление последней ячейки
def calc_last_cell(cell_prev: Cell, cell_current: Cell):
    x_cur = cell_current.x
    u_prev = cell_prev.pur_midl.u
    u_cur = cell_current.pur.u 
    p_prev = cell_prev.pur_midl.p
    r_prev = cell_prev.pur_midl.r  
    tau = layers[-1].tau

    def c(p, r):
        return math.sqrt(kappa * p / r)

    def p():
        return  p_prev * (1 - (kappa - 1) / (2 * c(p_prev, r_prev)) * (u_cur - u_prev)) ** (2 * kappa / (kappa - 1))

    def ro():
        return r_prev * (p() / p_prev) ** (1 / kappa)

    def x_n1():
        return x_cur + (u_cur + 2 * c(p(), ro()) / (kappa - 1)) * tau + 2 * c(p(), ro()) ** 2 / (kappa-1) * M / (S * p())* \
                                                            (1 - (1 + (kappa + 1) / (2 * c(p(), ro())) * p() * S / M * tau) ** (2 / (kappa + 1)))

    def u_n1():
        return u_cur + 2 * c(p(), ro()) / (kappa - 1) * (1 - (1 + (kappa + 1) / (2 * c(p(), ro())) * p() * S / M * tau) \
                                                                                        ** (-(kappa - 1) / (kappa + 1)))
    
    cell_current.x1 = x_n1()
    
    pur = PUR()
    pur.u = u_n1()

    p_prev_n1 = cell_prev.pur_midl_next.p
    u_cur_n1 = pur.u
    u_prev_n1 = cell_prev.pur_midl_next.u
    def p_n1():
        return p_prev_n1 * (1 - (kappa - 1) / (2 * c(p_prev, r_prev)) * (u_cur_n1 - u_prev_n1)) ** (2 * kappa / (kappa - 1))
    pur.p = p_n1()

    r_prev_n1 = cell_prev.pur_midl_next.r
    p_cur_n1 = pur.p
    def r_n1():
        return  r_prev_n1 * (p_cur_n1 / p_prev_n1)
    pur_r = r_n1()

    cell_current.pur_next(pur)


def calc_next_midl_cell():
        

def main():
    init()

    n = 0
    while layers[-1].x < L:
        new_layer = Layer()
        #Находим большие величины P U R
        #Установка половины PUR для первой ячейки
        nullPUR = PUR()
        firstCell = layers[-1].cells[1]
        nullPUR.p = firstCell.pur_midl.p
        nullPUR.u = -firstCell.pur_midl.u
        nullPUR.r = firstCell.pur_midl.r
        layers[-1].cells[0].pur_midl = nullPUR

        #От 1 до N - 1
        for i in range(1, N):
            calc_PUR(layers[-1].cells[i-1], layers[-1].cells[i])
        
        #N-я ячейка
        if(n == 0):
            prelast_PUR = layers[-1].cells[N-1].pur_midl
            layers[-1].cells[N].pur = prelast_PUR 
        else:
             calc_last_cell(layers[-1].cells[N - 1], layers[-1].cells[N])

        new_layer.h = layers[-1].cells[N].x1 / N



        


main()