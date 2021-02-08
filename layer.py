class Layer:
    cells = []
    h = 0
    tau = 0
    x = 0

    def append_cell(self, cell):
        self.cells.append(cell)