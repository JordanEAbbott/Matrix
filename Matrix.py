import itertools

def Levi_Civita(indices):
    parity = 1
    
    for index in range(1, len(indices)):
        for lower_index in range(0, index):
            try:
                parity *= (indices[index] - indices[lower_index]) / (abs(indices[index] - indices[lower_index]))
            except ZeroDivisionError:
                return 0
    
    return int(parity)

class Vector():
    
    def __init__(self, values, vec_type):
        self.dim = len(values)
        self.values = values
        self.type = vec_type
        
    def dot(self, other):
        
        return sum([a * b for a, b in zip(self.values, other.values)])
    
    def __mul__(self, multiplier):
        
        return Vector([multiplier * v for v in self.values], self.type)
    
    def __truediv__(self, factor):
        
        return Vector([v / factor for v in self.values], self.type)
    
    def __add__(self, other):
        
        if self.type != other.type:
            print("Vectors must be of same type")
            return
        
        return Vector([v + w for v, w in zip(self.values, other.values)], self.type)
    
    def __sub__(self, other):
        
        if self.type != other.type:
            print("Vectors must be of same type")
            return
        
        return Vector([v - w for v, w in zip(self.values, other.values)], self.type)
    
    def __str__(self):
        vector_string = '['
        if self.type == 'row':
            return str(self.values)
        else:
            for component in self.values:
                vector_string = vector_string + str(component) + ']\n['
            return vector_string[:-2]
        
class Matrix():
    
    def __init__(self, values, x_dim=False, y_dim=False):    
        if not x_dim or not y_dim:
            print("Put in matrix dimensions after the value list")
            
        v = 0
        self.x_dim = x_dim
        self.y_dim = y_dim
        self.m_dict = {}
        
        for i in range(1, self.x_dim + 1):
            for j in range(1, self.y_dim + 1):
                self.m_dict[(i, j)] = values[v]
                v += 1
                
    def det(self):
        if self.x_dim != self.y_dim:
            print("Cannot find the determinant of a non-square matrix")
            return
        det = 0
        
        indices = list(range(1, self.x_dim + 1))
        permuted = list(itertools.permutations(indices))
        for index_list in permuted:
            value = 1
            for i in range(1, self.x_dim + 1):
                value *= self.m_dict[(i, index_list[i - 1])]
            sign = Levi_Civita(index_list)
            det += sign * value
        
        return det
        
    def transpose(self):
        mT_dict = {}
        pseudo_dim = self.x_dim if self.x_dim > self.y_dim else self.y_dim
        
        for i in range(1, pseudo_dim + 1):
            for j in range(1, pseudo_dim + 1):
                try:
                    mT_dict[(i, j)] = self.m_dict[(j, i)]
                except:
                    pass
  
        return Matrix(list(mT_dict.values()), self.y_dim, self.x_dim)
    
    def adjugate(self):
        cof_dict = {}
        
        for i in range(1, self.x_dim + 1):
            for j in range(1, self.x_dim + 1):
                minor_values = [v for k, v in self.m_dict.items() if (k[0] != i and k[1] != j)]
                minor = Matrix(minor_values, self.x_dim - 1, self.y_dim - 1).det()
                cof_dict[(i, j)] = (-1)**(i + j) * minor
         
        cof_matrix = Matrix(list(cof_dict.values()), self.x_dim, self.y_dim)
        adj_matrix = cof_matrix.transpose()
        
        return adj_matrix
                    
    def inverse(self):
        det = self.det()
        if det == 0:
            print("Matrix is not invertible - 0 value determinant")
            return
        
        adj_dict = self.adjugate().m_dict
        for key in adj_dict:
            adj_dict[key] *= 1 / det
        
        return Matrix(list(adj_dict.values()), self.x_dim, self.y_dim)
    
    def dot(self, other):
        if self.y_dim != other.x_dim:
            print("Dimension Error: Can't multiply matrices")
            return
        
        dot_dict = {}
        
        for i in range(1, self.x_dim + 1):
            for k in range(1, other.y_dim + 1):
                row = Vector([self.m_dict[(i, j)] for j in range(1, self.y_dim + 1)], 'row')
                col = Vector([other.m_dict[(j, k)] for j in range(1, other.x_dim + 1)], 'col')
                dot_dict[(i, k)] = row.dot(col)
        
        return Matrix(list(dot_dict.values()), self.x_dim, other.y_dim)
    
        
    def swap_rows(self, row_1, row_2):
        
        temp = self[row_2].values
        self[row_2] = self[row_1].values
        self[row_1] = temp
        
        return
    
    def row_ech(self):
        
        pivot_row = 1
        pivot_column = 1
        
        while pivot_row <= self.x_dim and pivot_column <= self.y_dim:
            
            max_val = 0
            pivot = int
            for potential in range(1, self.x_dim + 1):
                if abs(self[potential, pivot_column]) > max_val:
                    max_val = abs(self[potential, pivot_column])
                    pivot = potential
                    
            if self[pivot, pivot_column] == 0:
                pivot_column += 1
            else:
                self.swap_rows(pivot_row, pivot)
                for i in range(pivot_row + 1, self.x_dim + 1):
                    fact = self[i, pivot_column] / self[pivot_row, pivot_column]
                    self[i, pivot_column] = 0
                    for j in range(pivot_column + 1, self.y_dim + 1):
                        self[i, j] = self[i, j] - self[pivot_row, j] * fact
                pivot_row += 1
                pivot_column += 1

        return
    
    def red_ech(self):
        
        lead = 0
        
        for i in range(0, self.x_dim):
            if self.y_dim <= lead:
                return
            j = i
            while self[j + 1, lead + 1] == 0:
                j = j + 1
                if self.x_dim == j:
                    j = i
                    lead = lead + 1
                    if self.y_dim == lead:
                        return
            if j != i:
                self.swap_rows(j + 1, i + 1)
            self[i + 1] = (self[i + 1] / self[i + 1, lead + 1]).values
            for k in range(0, self.x_dim):
                if k != i:
                    self[k + 1] = (self[k + 1] - self[i + 1] * self[k + 1, lead + 1]).values
            lead = lead + 1
        
        return
    
    def remove_row(self, row):
        
        return Matrix([v for k, v in self.m_dict.items() if k[0] != row], self.x_dim - 1, self.y_dim)
    
    def remove_column(self, column):
        
        m_list = [v for k, v in self.m_dict.items()]
        j = 0
        values = []
        for i in range(0, len(self.m_dict)):
            if i != (column - 1 + j * self.y_dim):
                values.append(m_list[i])
            else:
                j += 1
                
        return Matrix(values, self.x_dim, self.y_dim - 1)
    
    def decompose(self):  
        reduced = Matrix([v for k, v in self.m_dict.items()], self.x_dim, self.y_dim)
        reduced.red_ech()
        
        rows_to_remove = []
        for i in range(1, self.x_dim + 1):
            for j in range(1, self.y_dim + 1):
                if reduced[i, j] != 0:
                    break
                if j == self.y_dim:
                    rows_to_remove.append(i)
        tmp_F = reduced.remove_row(rows_to_remove[0])
        for row in rows_to_remove[1:]:
            tmp_F = tmp_F.remove_row(row)
            
        columns_to_remove = []
        leading_ones = [0 for i in range(0, self.y_dim)]
        for i in range(1, self.x_dim + 1):
            for j in range(1, self.y_dim + 1):
                if reduced[i, j] == 1:
                    leading_ones[j - 1] += 1
                    break
        for i in range(0, len(leading_ones)):
            if leading_ones[i] == 0:
                columns_to_remove.append(i + 1)
        tmp_C = self.remove_column(columns_to_remove[0])
        for column in columns_to_remove[1:]:
            tmp_C = tmp_C.remove_column(column)
           
        return  tmp_C, tmp_F
    
    def __getitem__(self, *keys):        
        
        if type(keys[0]) == int:
            return Vector([v for k, v in self.m_dict.items() if k[0] == keys[0]], 'row')
        
        return self.m_dict[(keys[0][0], keys[0][1])]
    
    def __setitem__(self, key, value):
        
        if type(key) == int:
            if len(value) != self.y_dim:
                print("Value list must be the length of a row to replace")
                return
            for i in range(1, self.y_dim + 1):
                self.m_dict[key, i] = value[i - 1]
            return

        self.m_dict[key[0], key[1]] = value
        
        return
    
    def __delitem__(self, keys):      
        return
    
    def __str__(self):
        matrix_string = '['
        
        for i in range(1, self.x_dim + 1):
            row = [self.m_dict[(i, j)] for j in range(1, self.y_dim + 1)]
            for element in row:
                matrix_string = matrix_string + str(element) + ', '
            matrix_string = matrix_string[:-2] + ']\n['
            
        return matrix_string[:-2]