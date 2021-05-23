import math
import cell
from decimal import Decimal

EPS = Decimal(10e-20)
kappa = Decimal(1.4)

def c(p, ro):
    return Decimal.sqrt(kappa * p / ro)

def Newtone(x, y, f1, f2, f_11, f_12, f_21, f_22) -> Decimal:
    while True:

        x1, y1 = x, y
        a11 = f_11(x1, y1)
        a12 = f_12(x1, y1)
        a21 = f_21(x1, y1)
        a22 = f_22(x1, y1)
        f11 = f1(x1, y1)
        f22 = f2(x1, y1)
        jacob = a11*a22 - a21*a12
        x = x1 - (a22*f11-a12*f22) / jacob
        y = y1 - (-a21*f11+a11*f22) / jacob
        if (Decimal.__abs__(f1(x, y) - f2(x, y)) < EPS):
            return x, y

def formul_1_2(cell_prev, cell_current) -> Decimal: #Проверил
    u_l = cell_prev.pur_midl.u
    u_r = cell_current.pur_midl.u
    p_l = cell_prev.pur_midl.p
    p_r = cell_current.pur_midl.p
    ro_l = cell_prev.pur_midl.r
    ro_r = cell_current.pur_midl.r

    def f1(u, p) -> Decimal:
        return u - u_r + 2 / (kappa - 1) * (1 - (p / p_r) ** ((kappa - 1) / (2 * kappa))) * c(p_r, ro_r)

    def f2(u, p) -> Decimal:
        return u - u_l - 2 / (kappa - 1) * (1 - (p / p_l) ** ((kappa - 1) / (2 * kappa))) * c(p_l, ro_l)  

    def f11(u, p) -> Decimal:
        return 1    

    def f12(u, p) -> Decimal:
        return (-1 / (kappa * p)) * c(p_r, ro_r) * (p / p_r) ** ((kappa - 1) / (2 * kappa))  

    def f21(u, p) -> Decimal:
        return 1  

    def f22(u, p) -> Decimal:
        return (1 / (kappa * p)) * c(p_l, ro_l) * (p / p_l) ** ((kappa - 1) / (2 * kappa))   

    return Newtone((u_l + u_r) / 2, (p_l + p_r) / 2, f1, f2, f11, f12, f21, f22)

def formul_1_3(cell_prev, cell_current) -> Decimal: #Проверил
    u_l = cell_prev.pur_midl.u
    u_r = cell_current.pur_midl.u
    p_l = cell_prev.pur_midl.p
    p_r = cell_current.pur_midl.p
    ro_l = cell_prev.pur_midl.r
    ro_r = cell_current.pur_midl.r

    def f1(u, p) -> Decimal:
        return u - u_r + 2 / (kappa - 1) * (1 - (p / p_r) ** ((kappa - 1) / (2 * kappa))) * c(p_r, ro_r)

    def f2(u ,p) -> Decimal:
        return u - u_r - (p - p_r) / Decimal.sqrt(ro_r * ((kappa + 1) / 2 * p + (kappa - 1) / 2 * p_r ))

    def f11(u, p) -> Decimal:
        return 1    

    def f12(u, p):
        return (-1 / (kappa * p)) * c(p_r, ro_r) * (p / p_r) ** ((kappa - 1) / (2 * kappa))  

    def f21(u, p) -> Decimal:
        return 1
    
    def f22(u , p) -> Decimal:
        return -(p_r - p) * (kappa + 1) / (4 * Decimal.sqrt(ro_r * (p_r * ((kappa - 1) / 2) + p * (kappa + 1) / 2)) * (p_r * ((kappa - 1) / 2) + p * (kappa + 1) / 2)) - \
            1 / (Decimal.sqrt(ro_r * (p_r * ((kappa - 1) / 2) + p * (kappa + 1) / 2)))

    return Newtone((u_l + u_r) / 2, (p_l + p_r) / 2, f1, f2, f11, f12, f21, f22)

def formul_1_4(cell_prev, cell_current) -> Decimal: #Проверил
    u_l = cell_prev.pur_midl.u
    u_r = cell_current.pur_midl.u
    p_l = cell_prev.pur_midl.p
    p_r = cell_current.pur_midl.p
    ro_l = cell_prev.pur_midl.r
    ro_r = cell_current.pur_midl.r

    def f1(u, p) -> Decimal:
        return u - u_r + 2 / (kappa - 1) * (1 - (p / p_r)**((kappa - 1)/(2 * kappa))) * c(p_r, ro_r)
    
    def f2(u,p) -> Decimal:
        return u - u_l + (p - p_l) / Decimal.sqrt(ro_l * ((kappa + 1) / 2 * p + (kappa - 1) / 2 * p_l ))
    
    def f11(u, p) -> Decimal:
        return 1    

    def f12(u, p) -> Decimal:
        return (-1 / (kappa * p)) * c(p_r, ro_r) * (p / p_r) ** ((kappa - 1) / (2 * kappa))  

    def f21(u, p) -> Decimal:
        return 1
    
    def f22(u , p) -> Decimal:
        return -(-p_l + p) * (kappa + 1) / (4 * Decimal.sqrt(ro_l * (p_l * ((kappa - 1) / 2) + p * (kappa + 1) / 2)) * (p_l * ((kappa - 1) / 2) + p * (kappa + 1) / 2)) + \
            1 / (Decimal.sqrt(ro_l * (p_l * ((kappa - 1) / 2) + p * (kappa + 1) / 2)))  

    return Newtone((u_l + u_r) / 2, (p_l + p_r) / 2, f1, f2, f11, f12, f21, f22)

def formul_3_4(cell_prev, cell_current) -> Decimal: #Проверил
    u_l = cell_prev.pur_midl.u
    u_r = cell_current.pur_midl.u
    p_l = cell_prev.pur_midl.p
    p_r = cell_current.pur_midl.p
    ro_l = cell_prev.pur_midl.r
    ro_r = cell_current.pur_midl.r
    
    def f1(u ,p) -> Decimal:
        return u - u_r - (p - p_r) / Decimal.sqrt(ro_r * ((kappa + 1) / 2 * p + (kappa - 1) / 2 * p_r ))

    def f2(u,p) -> Decimal:
        return u - u_l + (p - p_l) / Decimal.sqrt(ro_l * ((kappa + 1) / 2 * p + (kappa - 1) / 2 * p_l ))

    def f11(u, p) -> Decimal:
        return 1
    
    def f12(u , p) -> Decimal:
        return -(p_r - p) * (kappa + 1) / (4 * Decimal.sqrt(ro_r * (p_r * ((kappa - 1) / 2) + p * (kappa + 1) / 2)) * (p_r * ((kappa - 1) / 2) + p * (kappa + 1) / 2)) - \
            1 / (Decimal.sqrt(ro_r * (p_r * ((kappa - 1) / 2) + p * (kappa + 1) / 2)))

    def f21(u, p) -> Decimal:
        return 1
    
    def f22(u , p) -> Decimal:
        return -(-p_l + p) * (kappa + 1) / (4 * Decimal.sqrt(ro_l * (p_l * ((kappa - 1) / 2) + p * (kappa + 1) / 2)) * (p_l * ((kappa - 1) / 2) + p * (kappa + 1) / 2)) + \
            1 / (Decimal.sqrt(ro_l * (p_l * ((kappa - 1) / 2) + p * (kappa + 1) / 2))) 

    return Newtone((u_l + u_r) / 2, (p_l + p_r) / 2, f1, f2, f11, f12, f21, f22)