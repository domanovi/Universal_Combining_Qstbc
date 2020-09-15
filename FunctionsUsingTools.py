from Tools import *
import pandas as pd

from sympy import *
import math
import copy


# Equivalent MISO functions
# ===================================================================================================


def produce_abba_qstbc(indices):
    """calculates a matrix of information variables indexed by 'indices'
     forming an ABBA-QSTBC structure"""
    len_indices = len(indices)
    if not math.log2(len_indices).is_integer() or len_indices < 2:
        return
    if len_indices == 2:
        return Block([[Element("x_{" + str(indices[0]) + "}"), Element("x_{" + str(indices[1]) + "}")],
                      [Element("x_{" + str(indices[1]) + "}*"), Element("-x_{" + str(indices[0]) + "}*")]])
    else:
        a = produce_abba_qstbc(indices[0: int(len_indices / 2)])
        b = produce_abba_qstbc(indices[int(len_indices / 2):])
        return Block([[a, b], [b, a]])


def produce_abba_qstbc_of_order(n):
    """calculates a matrix of n information variables
     forming an ABBA-QSTBC structure"""
    return produce_abba_qstbc(range(n))


def produce_ordered_equivalent_channel(n):
    """produces the channel coefficient vector in the MISO equivalent channel"""
    h = Block([[Element(chr(0x0125) + "_{" + str(ind) + "}")] for ind in range(n)])
    return h


def produce_equivalent_channel_vec(p, q, rec_indices=None):
    """Produces the the equivalent channel vector in terms of the actual channel coefficients
    h_{i, j} where 'i' is the index of the transmitter and 'j' is the index of the receiever."""
    tx_num, rx_num = 2 ** p, 2 ** q
    # ------ 2 parameters call ------
    if rec_indices is None:
        return produce_equivalent_channel_vec(p, q, range(rx_num))
    # ------ 3 parameters call ------
    else:
        if len(rec_indices) != rx_num:
            raise ValueError("length of 'rec_indices' vector must be equal to 2^q")
        # For 2 receivers
        if q == 1:
            elem_list = []
            for tx_ind in range(tx_num):
                elem_list += [[Element("h_{" + str(tx_ind) + "-" + str(rec_indices[0]) + "}")]]
                elem_list += [[Element("h_{" + str(tx_ind) + "-" + str(rec_indices[1]) + "}*")]]
            return Block(elem_list)
        # For more than 2 receivers
        else:
            h_hat_low = produce_equivalent_channel_vec(p, q - 1, rec_indices[0: int(rx_num / 2)])
            h_hat_high = produce_equivalent_channel_vec(p, q - 1, rec_indices[int(rx_num / 2):])
            return Block([[h_hat_low], [h_hat_high]])


# Actual MIMO channel
# ===================================================================================================


def produce_transmission(p, q, indices=None):
    """calculates a transmission matrix of information variables indexed by 'indices' for 2^p transmit
    and 2^q receive antennas MIMO scheme, which leads to ABBA-QSTBC equivalent MISO scheme.
    The length of 'indices' must be 2^(p+q).
    If 'indices' is not specified it will be set to be the vector of indices from 0 to 2^(p+q)-1."""
    # ------ 2 parameters call -------
    if indices is None:
        return produce_transmission(p, q, range(2 ** (p + q)))
    # ------ 3 parameters call -------
    len_indices = len(indices)
    if len_indices != 2 ** (p + q):
        raise ValueError("length of 'indices' vector must be equal to 2 ^ (p+q)")
    # For 2 receivers
    if q == 1:
        if p == 0:
            return Block([[Element("x_{" + str(indices[0]) + "}")],
                          [Element("x_{" + str(indices[1]) + "}*")]])
        else:
            c_low = produce_transmission(p - 1, q, indices[0: int(len_indices / 2)])
            c_high = produce_transmission(p - 1, q, indices[int(len_indices / 2):])
            return Block([[c_low, c_high], [c_high, c_low]])
    # For more than 2 receivers
    else:
        c_low = produce_transmission(p, q - 1, indices[0: int(len_indices / 2)])
        c_high = produce_transmission(p, q - 1, indices[int(len_indices / 2):])
        return Block([[c_low], [c_high]])


def produce_operation_matrix(p, q):
    """calculates a operation matrix for 2^p transmitters and 2^q  receivers MIMO scheme
     which leads to ABBA-QSTBC equivalent MISO scheme"""
    tx_num, rx_num = 2 ** p, 2 ** q
    # For 2 receivers
    if q == 1:
        if p == 0:
            return Matrix([[Symbol("+"), Symbol("-*")],
                           [Symbol("+"), Symbol("*")]])
        else:
            op_mat = produce_operation_matrix(p - 1, q)
            return Matrix([[op_mat], [op_mat]])
    # For more than 2 receivers
    else:
        op_mat = produce_operation_matrix(p, q - 1)
        return Matrix([[op_mat, op_mat], [op_mat, op_mat]])


def build_channel_matrix(p, q):
    """Builds a matrix of actual channel coefficients for qxp MIMO scheme, such that in row 'i' and colomn
    'j' the symbol 'h_{i,j}' will be set."""
    tx_num, rx_num = 2 ** p, 2 ** q
    return Block([[Element("h_{" + str(i) + "-" + str(j) + "}") for j in range(rx_num)] for i in range(tx_num)])


def operate_matrix(op, mat):
    """Returns the matrix 'mat' after operating the operations specified in 'op' matrix on each element."""
    if op.shape != mat.shape:
        raise ValueError('operation matrix must be at the same shape of the matrix itself')
    result = copy.deepcopy(mat)
    [row_num, col_num] = mat.shape
    for i in range(row_num):
        for j in range(col_num):
            if str(op[i, j]) == "+":
                result[i, j] = mat[i, j]
            elif str(op[i, j]) == "-":
                result[i, j] = -mat[i, j]
            elif str(op[i, j]) == "*":
                result[i, j] = mat[i, j].conj()
            elif str(op[i, j]) == "-*":
                result[i, j] = mat[i, j].conj().minus()
            else:
                raise ValueError('operations must be one of {-,+,*,=*}')
    return result


def cross_addition(mat):
    """Calculates the 'Cross Addition' operation on the matrix.
    The number of rows must be even, and the number of colomn must be a power of 2."""
    row_num, col_num = mat.shape
    q = math.log2(col_num)
    if row_num % 2 != 0 or not q.is_integer():
        raise ValueError('matrix dimensions are not legit: must be with even number of rows, and power of'
                         '2 number of colomns')
    # For 2 receivers
    if q == 1:
        result = Block([[Element("0")]]*row_num)
        for i in range(row_num):
            if i % 2 == 0:
                result[i, 0] = mat[i, 0] + mat[i + 1, 1]
            else:
                result[i, 0] = mat[i, 0] + mat[i - 1, 1]
        return result
    # For more than 2 receivers
    else:
        a = mat[0: int(row_num / 2), 0: int(col_num / 2)]
        b = mat[0: int(row_num / 2), int(col_num / 2):]
        c = mat[int(row_num / 2):, 0:int(col_num / 2)]
        d = mat[int(row_num / 2):, int(col_num / 2):]
        return Block([[cross_addition(a) + cross_addition(d)],
                      [cross_addition(c) + cross_addition(b)]])


# Constant Matrices Representation
# ===================================================================================================
def transmission_column_representation(trans_matrix):
    """ The function gets a transmission matrix and returns a list with Column Representation of it.
    Each column of the 'trans_matrix' can be represented as a multiplication of the corresponding
    Column Representation matrix of constants and the vectors composed of the information symbols and
    their conjugates, for example : [x_1, x_2, x_3, x_1*, x_2*, x_3*] ^ T"""
    times, tx_num = trans_matrix.shape
    # find max index
    max_index = max([max([int(elem.symbol.indices_list[0])
                          for elem in trans_matrix[i]]) for i in range(times)])
    # Build matrices list
    col_num = 2 * (max_index + 1)
    result = []
    for tx_ind in range(tx_num):
        tmp_result = np.zeros((times, col_num))
        for time_ind in range(times):
            index = int(trans_matrix[time_ind, tx_ind].symbol.indices_list[0])
            elem_conj = trans_matrix[time_ind, tx_ind].conjugate
            elem_sign = trans_matrix[time_ind, tx_ind].sign
            if elem_conj == "*":
                index += max_index + 1
            # update the relevant matrix entry
            tmp_result[time_ind, index] = 1 * (elem_sign == '+') + (-1) * (elem_sign == '-')
        result.append(tmp_result)
    return result


def combining_matrix_representation(comb_matrix):
    rows_num = comb_matrix.shape[0]
    # find maximum receiver and time index
    max_time = max(
        [max([int(elem.symbol.indices_list[TIME_INDEX]) for elem in comb_matrix[row_ind, 0].block[0]]) for row_ind in
         range(rows_num)])
    max_rec = max(
        [max([int(elem.symbol.indices_list[REC_INDEX]) for elem in comb_matrix[row_ind, 0].block[0]]) for row_ind in
         range(rows_num)])
    # Build matrices list
    col_num = 2 * (max_time + 1)
    result = []
    for rx_ind in range(max_rec + 1):
        tmp_result = np.zeros((rows_num, col_num))
        # for each row find an element with 'rx_ind'
        for row_ind in range(rows_num):
            for elem in comb_matrix[row_ind, 0].block[0]:
                if int(elem.symbol.indices_list[REC_INDEX]) == rx_ind:
                    elem_time = int(elem.symbol.indices_list[TIME_INDEX])
                    elem_conj = elem.conjugate
                    elem_sign = elem.sign
                    if elem_conj == '*':
                        time_index = elem_time + max_time + 1
                    else:
                        time_index = elem_time
                    tmp_result[row_ind, time_index] = 1 * (elem_sign == '+') + (-1) * (elem_sign == '-')
        result.append(tmp_result)
    return result


# Only for Presentation
# ===================================================================================================
TIME_INDEX = 0
REC_INDEX = 1


class TimedElement(Element):

    def __init__(self, index, time):
        self.index = index
        self.time = time
        Element.__init__(self, elem_string="r_{" + str(time) + "-" + str(index) + "}")

    def __repr__(self):
        sign_str = self.sign
        conj_str = chr(SUPERSCRIPT_STAR) * (self.conjugate == '*') +\
                   '' * (self.conjugate != '*')
        if self.index == 'none':
            rec_num_str = 'none'
        else:
            rec_num_str = ''
            for char in str(self.index):  # first index
                rec_num_str += chr(SUBSCRIPT_0 + int(char))
        return 'r' + rec_num_str + conj_str + '(' + str(self.time) + ')'


def produce_combining_matrix(p, q):
    op_mat = produce_operation_matrix(p, q)
    num_rec = 2 ** q
    times = (2 ** q) * (2 ** p)
    comb_mat = Block([[TimedElement(0, 0)] * num_rec] * times)
    for time_ind in range(times):
        for rec_ind in range(num_rec):
            comb_mat[time_ind, rec_ind] = TimedElement(rec_ind, time_ind)
    return cross_addition(operate_matrix(op_mat, comb_mat))


# Printing
# ===================================================================================================
class SmartStr(str):
    """ The class SmartStr is a string implementing a special '__add__' dunction
    that smartly adds multiple line strings horizontally"""
    def __add__(self, other):
        """ The function gets 2 strings, each has multiple lines saperated by '\n'.
        It returns a single string which represents an horizontal concatenation of the 2 strings
        where each line includes the corresponding lines at both str1 and str2."""
        list1, list2 = str(self).split('\n'), str(other).split('\n')
        # Append lines to the shorter list so the number of lines would be equal
        if len(list1) < len(list2):
            list1 += [''] * (len(list2) - len(list1))
        else:
            list2 += [''] * (len(list1) - len(list2))
        # fill each string in each list with spaces so each line would be with the same length
        max_line1, max_line2 = max([len(line) for line in list1]), max([len(line) for line in list2])
        for ind, line in enumerate(list1):
            list1[ind] += ' ' * (max_line1 - len(line))
        for ind, line in enumerate(list2):
            list2[ind] += ' ' * (max_line2 - len(line))
        # Concatenate the lists line-wise
        result = ''
        for line_ind in range(len(list1)):
            result += list1[line_ind] + list2[line_ind]
            if line_ind != len(list1) - 1:
                result += '\n'
        return SmartStr(result)


def print_information_symbol_vector_str(max_ind):
    return str(Block("//".join(["x_{{{0}}}".format(ind) for ind in range(max_ind)]) + "//"
                 + "//".join(["x_{{{0}}}*".format(ind) for ind in range(max_ind)])).transpose()) + chr(SUPERSCRIPT_T)


SUBSCRIPT_J = 0x2C7C
SUPERSCRIPT_T = 0x1D40

def print_received_stream_rx_i_str(max_time):
    return "[[ "  + "   ".join(["r{}({})".format(chr(SUBSCRIPT_J), time) for time in range(max_time)]) + "   "\
           + "   ".join(["r{}({}){}".format(chr(SUBSCRIPT_J), time, chr(SUPERSCRIPT_STAR)) for time in range(max_time)])\
           + " ]]" + chr(SUPERSCRIPT_T)


def present_qstbc_scheme(p, q):
    n_t, n_r = 2 ** p, 2 ** q
    print(" ================================= Target ==================================")
    print("Nt = " + str(n_t))
    print("Nr = " + str(n_r))
    print("Target QSTBC:")
    stbc = produce_abba_qstbc_of_order(2 ** p * 2 ** q)
    print(stbc)
    print("Desired Stream:")
    h_hat = produce_ordered_equivalent_channel(n_t * n_r)
    s_hat = np.matmul(stbc, h_hat)
    print(SmartStr(chr(0x015D) + " = ") + SmartStr(stbc) + SmartStr(h_hat) + "  =  " + SmartStr(s_hat))
    print("\n ================================= Result ==================================")
    print("Transmission Matrix:")
    trans_matrix = produce_transmission(p, q)
    print(trans_matrix)
    print("Combining Scheme:")
    comb_scheme = produce_combining_matrix(p, q)
    print(comb_scheme)
    print("Channel coefficients transformation:")
    h_transform = produce_equivalent_channel_vec(p, q)
    print(SmartStr(h_hat) + " = " + SmartStr(h_transform))
    print("\n =============================== Check Result ==============================")
    print("Received Signals Matrix:")
    H = build_channel_matrix(p, q)
    r = np.matmul(trans_matrix, H)
    print(SmartStr("\nR    =    ") + SmartStr(pandas.DataFrame(r.view(np.ndarray)).to_string()))
    print("Received Stream (after combining):")
    op_mat = produce_operation_matrix(p, q)
    s = cross_addition(operate_matrix(op_mat, r))
    print(SmartStr("s = ") + SmartStr(s))
    print("Compare 'Desired Stream' and 'Received Stream':")
    s_hat_to_compare = np.matmul(stbc, h_transform)
    result = (s == s_hat_to_compare)
    print(result)
    if any(result):
        print("\n ========================== Comparison Succeeded !!! ==========================")
    else:
        print("\n ========================== Comparison Failed ===========================")


def present_const_matrices(p, q):
    n_t, n_r = 2 ** p, 2 ** q
    print("Nt = " + str(n_t))
    print("Nr = " + str(n_r))
    print("\n =============================== Transmission Matrices ==============================")
    trans_matrix = produce_transmission(p, q)
    trans_matrices = transmission_column_representation(trans_matrix)
    print("Each matrix is multiplied from right by the information symbols vector:")
    print(print_information_symbol_vector_str(int(trans_matrices[0].shape[1]/2)))
    trans_str = SmartStr('')
    for ind, mat in enumerate(trans_matrices):
        trans_str += SmartStr('Tx ' + str(ind) + ':\n' + str(mat) + '    ')
    print(trans_str)
    print("\n =============================== Combining Matrices ==============================")
    comb_matrix = produce_combining_matrix(p, q)
    comb_matrices = combining_matrix_representation(comb_matrix)
    print("The stream received in rx i (and its conjugates terms):")
    print(print_received_stream_rx_i_str(int(comb_matrices[0].shape[1]/2)))
    print("The stream received in rx i (and its conjugates terms) is multiplied\n"
          " from left by the corresponding combining matrix and then the results are summed up.")
    comb_str = SmartStr('')
    for ind, mat in enumerate(comb_matrices):
        comb_str += SmartStr('Rx ' + str(ind) + ':\n' + str(mat) + '    ')
    print(comb_str)