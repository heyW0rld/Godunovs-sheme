from decimal import Decimal

from pur import PUR


class Cell:
    x = Decimal(0)

    def __init__(self):
        self.pur = PUR()
        self.pur_midl = PUR()
   
    # def set_values(self, p, ro, u):
    #     self.pur.p = p
    #     self.pur.ro = ro
    #     self.pur.u = u

    # def set_midl_values(self, p, ro, u):
    #     self.pur_midl.p = p
    #     self.pur_midl.ro = ro
    #     self.pur_midl.u = u

    # def set_next_values(self, p, ro, u):
    #     self.pur_next.p = p
    #     self.rur_next.ro = ro
    #     self.pur_next.u = u
    
    # def set_next_midl_values(self, p, ro, u):
    #     self.pur_midl_next.p = p
    #     self.pur_midl_next.ro = ro
    #     self.pur_midl_next.u = u

