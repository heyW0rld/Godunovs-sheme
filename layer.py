from decimal import Decimal
from cell import Cell

class Layer:
    h = Decimal(0)
    tau = Decimal(0)
    x = Decimal(0)

    def __init__(self):
        self.cells = []

    def append_cell(self, cell):
        self.cells.append(cell)