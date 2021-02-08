import math
import cell

EPS = 10e-3
kappa = 1.4

def c(p, ro):
    return math.sqrt(kappa * p / ro)

def Newtone(x, y, f1, f2, f_11, f_12, f_21, f_22):
    while True:
        x1, y1 = x, y
        a11 = f_11(x1, y1)
        a12 = f_12(x1, y1)
        a21 = f_21(x1, y1)
        a22 = f_22(x1, y1)
        f11 = f1(x1, y1)
        f22 = f2(x1, y1)
        jacob = a11*a22 - a21*a12
        x = x1 - (a22*f11-a12*f22)/jacob
        y = y1 - (-a21*f11+a11*f22) / jacob
        if (abs(f1(x, y) - f2(x, y)) < EPS):
            return x, y

def formul_1_2(cell_prev, cell_current):
    u_l = cell.u_prev
    u_r = cell.u_next
    p_l = cell.p_prev
    p_r = cell.p_next
    ro_l = cell.ro_prev
    ro_r = cell.ro_next

    def f1(u, p):
       return u - u_r + 2/(kappa - 1) * (1 - (p/p_r)**((kappa - 1)/(2*kappa)))*c(p_r, ro_r)

    def f2(u, p):
        return u - u_l - 2/(kappa - 1) * (1 - (p/p_l)**((kappa - 1)/(2*kappa)))*c(p_l, ro_l)  

    def f11(u, p):
        return 1    

    def f12(u, p):
        return -1/kappa * c(p_r, ro_r) * p**((kappa-1)/(2*kappa) - 1) / p_r**((kappa-1)/(2*kappa)) 

    def f21(u, p):
        return 1  

    def f22(u, p):
        return 1/kappa * c(p_r, ro_r) * p**((kappa-1)/(2*kappa) - 1) / p_r**((kappa-1)/(2*kappa))   

    return Newtone(cell.u_next, cell.p_next, f1, f2, f11, f12, f21, f22)

def formul_1_3(cell_prev, cell_current):
    u_l = cell.u_prev
    u_r = cell.u_next
    p_l = cell.p_prev
    p_r = cell.p_next
    ro_l = cell.ro_prev
    ro_r = cell.ro_next

    def f1(u, p):
       return u - u_r + 2/(kappa - 1) * (1 - (p/p_r)**((kappa - 1)/(2*kappa)))*c(p_r, ro_r)

    def f2(u ,p):
        return u - u_r - (p - p_r)/math.sqrt(p_r * ((kappa + 1) / 2 * p + (kappa - 1) / 2 * p_r ))

    def f11(u, p):
        return 1    

    def f12(u, p):
        return -1/kappa * c(p_r, ro_r) * p**((kappa-1)/(2*kappa) - 1) / p_r**((kappa-1)/(2*kappa)) 

    def f21(u, p):
        return 1
    
    def f22(u , p):
        return (2 * p_r**2*(kappa-1)+(kappa+1)*p_r*p) / (math.sqrt(2) * (p_r**2*(kappa - 1) + p_r * (kappa + 1) * p))**(3/2)

    return Newtone(cell.u_next, cell.p_next, f1, f2, f11, f12, f21, f22)

def formul_1_4(cell_prev, cell_current):
    u_l = cell.u_prev
    u_r = cell.u_next
    p_l = cell.p_prev
    p_r = cell.p_next
    ro_l = cell.ro_prev
    ro_r = cell.ro_next

    def f1(u, p):
        return u - u_r + 2/(kappa - 1) * (1 - (p/p_r)**((kappa - 1)/(2*kappa)))*c(p_r, ro_r)
    
    def f2(u,p):
        return u - u_l + (p - p_l)/math.sqrt(p_l * ((kappa + 1) / 2 * p + (kappa - 1) / 2 * p_l ))
    
    def f11(u, p):
        return 1    

    def f12(u, p):
        return -1/kappa * c(p_r, ro_r) * p**((kappa-1)/(2*kappa) - 1) / p_r**((kappa-1)/(2*kappa)) 

    def f21(u, p):
        return 1
    
    def f22(u , p):
        return (2 * p_l**2*(kappa-1)+(kappa+1)*p_l*p) / (math.sqrt(2) * (p_l**2*(kappa - 1) + p_l * (kappa + 1) * p))**(3/2)

    return Newtone(cell.u_next, cell.p_next, f1, f2, f11, f12, f21, f22)

def formul_3_4(cell_prev, cell_current):
    u_l = cell.u_prev
    u_r = cell.u_next
    p_l = cell.p_prev
    p_r = cell.p_next
    ro_l = cell.ro_prev
    ro_r = cell.ro_next
    
    def f1(u ,p):
        return u - u_r - (p - p_r)/math.sqrt(p_r * ((kappa + 1) / 2 * p + (kappa - 1) / 2 * p_r ))

    def f2(u,p):
        return u - u_l + (p - p_l)/math.sqrt(p_l * ((kappa + 1) / 2 * p + (kappa - 1) / 2 * p_l ))

    def f11(u, p):
        return 1
    
    def f12(u , p):
        return (2 * p_r**2*(kappa-1)+(kappa+1)*p_r*p) / (math.sqrt(2) * (p_r**2*(kappa - 1) + p_r * (kappa + 1) * p))**(3/2)

    def f21(u, p):
        return 1
    
    def f22(u , p):
        return (2 * p_l**2*(kappa-1)+(kappa+1)*p_l*p) / (math.sqrt(2) * (p_l**2*(kappa - 1) + p_l * (kappa + 1) * p))**(3/2)

    return Newtone(cell.u_next, cell.p_next, f1, f2, f11, f12, f21, f22)