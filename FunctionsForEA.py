from Tools import *
import pandas as pd

from sympy import *
import math
import copy


# Equivalent MISO functions
# ===================================================================================================


def produce_generalized_ea_qstbc(indices):
    """calculates a matrix of information variables indexed by 'indices'
     forming an extended EA-QSTBC structure.
    'indices' must be with power of 2 number of elements."""
    len_indices = len(indices)
    if not math.log2(len_indices).is_integer() or len_indices < 1:
        return
    if len_indices == 1:
        return Block([[Element("x_{" + str(indices[0]) + "}")]])
    else:
        a = produce_generalized_ea_qstbc(indices[0: int(len_indices / 2)])
        b = produce_generalized_ea_qstbc(indices[int(len_indices / 2):])
        return Block([[a, b], [b.conj().minus(), a.conj()]])


def produce_ea_qstbc_of_order(n):
    """calculates a matrix of n information variables
     forming an n x n EQ-QSTBC structure.
     n must be integer power of 2."""
    return produce_generalized_ea_qstbc(range(n))


def produce_punctured_ea_qstbc(n_tild, n, punctured_columns):
    """calculates a matrix for n transmit antennas of n_tild information variables,
     where n_tild is the closest power of 2 which is greater than n.
     The result matrix is built such that the n_tild x n_tild EA-QSTBC is built first,
     and then the columns indexed by 'punctured_columns' are punctured."""

    # check input
    if len(punctured_columns) != n_tild - n:
        raise ValueError("the punctured_columns vector must be of size 'n_tild - n'")
    if not all(isinstance(index, int) and 0 <= index < n_tild for index in punctured_columns):
        raise ValueError("Some the punctured indices are not in the right range")
    # produce matrix
    left_columns = list(range(n_tild))
    for index in punctured_columns:
        left_columns.remove(index)
    ea_matrix = produce_ea_qstbc_of_order(n_tild)
    return ea_matrix[:, left_columns]


# ------------------------------------------------------------------------------


def produce_equivalent_channel_vector(n):
    """produces the channel coefficient vector in the MISO equivalent channel"""
    h = Block([[Element(chr(0x0125) + "_{" + str(ind) + "}")] for ind in range(n)])
    return h


def produce_extended_reordering_transformation(p, q, rec_indices=None):
    """Produces the the equivalent channel vector in terms of the actual channel coefficients
    h_{i, j} where 'i' is the index of the transmitter and 'j' is the index of the receiver."""
    tx_num, rx_num = 2 ** p, 2 ** q
    # ------ 2 parameters call ------
    if rec_indices is None:
        return produce_extended_reordering_transformation(p, q, range(rx_num))
    # ------ 3 parameters call ------
    else:
        if len(rec_indices) != rx_num:
            raise ValueError("length of 'rec_indices' vector must be equal to 2^q")
        # For 1 receiver
        if q == 0:
            elem_list = []
            for tx_ind in range(tx_num):
                elem_list += [[Element("h_{" + str(tx_ind) + "-" + str(rec_indices[0]) + "}")]]
            return Block(elem_list)
        # For more than 2 receivers
        else:
            h_hat_low = produce_extended_reordering_transformation(p, q - 1, rec_indices[0: int(rx_num / 2)])
            h_hat_high = produce_extended_reordering_transformation(p, q - 1, rec_indices[int(rx_num / 2):]).conj()
            return Block([[h_hat_low], [h_hat_high]])


def produce_punctured_reordering_transformation(n_t, n_r, punctured_tx, punctured_rx):
    """The function produces the equivalent channel vector such that all the coefficients that involve
    indices from 'punctured_tx' or 'punctured_rx' are removed"""
    p, q = math.ceil(math.log2(n_t)), math.ceil(math.log2(n_r))
    h_eq = produce_extended_reordering_transformation(p, q)
    # h_eq_punctured = h_eq[h_eq]
    punctured_indices = []
    for raw_ind, raw in enumerate(h_eq):
        if int(raw[0].symbol.indices_list[0]) in punctured_tx or int(raw[0].symbol.indices_list[1]) in punctured_rx:
            punctured_indices.append(raw_ind)
    h_eq_punctured = np.delete(h_eq, punctured_indices, 0)
    return h_eq_punctured, punctured_indices


def produce_reordering_transformation(n_t, n_r):
    """Produces the reordering transformation such that the last tx or rx antennas are punctured"""
    n_t_tild, n_r_tild = 2 ** math.ceil(math.log2(n_t)), 2 ** math.ceil(math.log2(n_r))
    return produce_punctured_reordering_transformation(n_t, n_r,
                                                       range(n_t, n_t_tild), range(n_r, n_r_tild))


# MIMO channel
# ===================================================================================================


def produce_extended_transmission(p, q, indices=None):
    """calculates a transmission matrix of information variables indexed by 'indices' for 2^p transmit
    and 2^q receive antennas MIMO scheme, which leads to EA-QSTBC equivalent MISO scheme.
    The length of 'indices' must be 2^(p+q).
    If 'indices' is not specified it will be set to be the vector of indices from 0 to 2^(p+q)-1."""
    # ------ 2 parameters call -------
    if indices is None:
        return produce_extended_transmission(p, q, range(2 ** (p + q)))
    # ------ 3 parameters call -------
    len_indices = len(indices)
    if len_indices != 2 ** (p + q):
        raise ValueError("length of 'indices' vector must be equal to 2 ^ (p+q)")
    # For 1 receivers
    if q == 0:
        return produce_generalized_ea_qstbc(indices)
    # For more than 1 receiver
    else:
        c_low = produce_extended_transmission(p, q - 1, indices[0: int(len_indices / 2)])
        c_high = produce_extended_transmission(p, q - 1, indices[int(len_indices / 2):]).conj()
        return Block([[c_low], [c_high]])


def produce_punctured_transmission(n_t, n_r, punctured_tx):
    """ Builds a transmission scheme for general n_t x n_r MIMO channel"""
    p, q = math.ceil(math.log2(n_t)), math.ceil(math.log2(n_r))
    trans_matrix = produce_extended_transmission(p, q)
    trans_matrix_punctured = np.delete(trans_matrix, punctured_tx, 1)
    return trans_matrix_punctured


def produce_transmission(n_t, n_r):
    """ Builds the general transmission while puncturing the last tx antennas """
    n_t_tild = 2 ** math.ceil(math.log2(n_t))
    return produce_punctured_transmission(n_t, n_r, range(n_t, n_t_tild))


# ------------------------------------------------------------------------------

def build_channel_matrix(n_t, n_r):
    """Builds a matrix of actual channel coefficients for n_t x n_r MIMO scheme, such that in row 'i' and colomn
    'j' the symbol 'h_{i,j}' will be set."""
    return Block([[Element("h_{" + str(i) + "-" + str(j) + "}") for j in range(n_r)] for i in range(n_t)])


def do_extended_combining(mat):
    """Calculates the combined vector after operating the combining scheme on the matrix.
    The number of rows must be even, and the number of colomn must be a power of 2."""
    row_num, col_num = mat.shape
    q = math.log2(col_num)
    if (row_num > 1 and row_num % 2 != 0) or not q.is_integer():
        raise ValueError('matrix dimensions are not legit: must be with even number of rows (or exactly 1), '
                         'and power of 2 number of colomns')
    # For 1 receiver
    if q == 0:
        sums_list = []
        for row_ind in range(row_num):
            elem = mat[row_ind, 0]
            if not isinstance(elem, Element):
                sums_list += [[elem]]
            else:
                sums_list += [[mSum(elem)]]
        return Block(sums_list)
    # For more than 1 receivers
    else:
        a = mat[0: int(row_num / 2), 0: int(col_num / 2)]
        b = mat[0: int(row_num / 2), int(col_num / 2):].conj()
        c = mat[int(row_num / 2):, 0:int(col_num / 2)].minus()
        d = mat[int(row_num / 2):, int(col_num / 2):].conj()
        return Block([[do_extended_combining(a) + do_extended_combining(d)],
                      [do_extended_combining(c) + do_extended_combining(b)]])


def do_punctured_combining(mat, punctured_rx):
    """The function combine the received matrix 'mat' to generate a combined vector, while using the puncturing
    over the receivers with 'punctured_rx' indices."""
    typ = type(mat[0, 0])
    times, _ = mat.shape
    zeros_col = Block([[typ(is_zero=True)]] * times)
    for index in punctured_rx:
        mat = np.hstack((mat[:, :index], zeros_col, mat[:, index:])).view(Block)
    return do_extended_combining(mat)


def do_combining(mat):
    """ Operates the combining scheme on a received matrix induced from a general MIMO system """
    times_num, n_r = mat.shape
    n_r_tild = 2 ** math.ceil(math.log2(n_r))
    return do_punctured_combining(mat, range(n_r, n_r_tild))


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

    def __init__(self, index=-1, time=-1, is_zero=False):
        self.is_zero = is_zero
        self.index = index
        self.time = time
        Element.__init__(self, elem_string="r_{" + str(time) + "-" + str(index) + "}", is_zero=is_zero)

    def __repr__(self):
        if self.is_zero is True:
            return ' 0 '
        sign_str = self.sign
        sign_str = '' * (self.sign == '+') + \
                   '-' * (self.conjugate == '-')
        conj_str = chr(SUPERSCRIPT_STAR) * (self.conjugate == '*') + \
                   '' * (self.conjugate != '*')
        if self.index == 'none':
            rec_num_str = 'none'
        else:
            rec_num_str = ''
            for char in str(self.index):  # first index
                rec_num_str += chr(SUBSCRIPT_0 + int(char))
        return sign_str + 'r' + rec_num_str + conj_str + '(' + str(self.time) + ')'


def produce_combining_matrix(n_t, n_r):
    n_t_tild, n_r_tild = 2 ** math.ceil(math.log2(n_t)), 2 ** math.ceil(math.log2(n_r))
    times = n_t_tild * n_r_tild
    comb_mat = Block([[TimedElement(0, 0)] * n_r] * times)
    for time_ind in range(times):
        for rec_ind in range(n_r):
            comb_mat[time_ind, rec_ind] = TimedElement(rec_ind, time_ind)
    return do_combining(comb_mat)


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
                     + "//".join(["x_{{{0}}}*".format(ind) for ind in range(max_ind)])).transpose()) + chr(
        SUPERSCRIPT_T)


SUBSCRIPT_J = 0x2C7C
SUPERSCRIPT_T = 0x1D40


def print_received_stream_rx_i_str(max_time):
    return "[[ " + "   ".join(["r{}({})".format(chr(SUBSCRIPT_J), time) for time in range(max_time)]) + "   " \
           + "   ".join(["r{}({}){}".format(chr(SUBSCRIPT_J), time, chr(SUPERSCRIPT_STAR)) for time in range(max_time)]) \
           + " ]]" + chr(SUPERSCRIPT_T)


def present_qstbc_scheme(n_t, n_r):
    if n_t <= 0 or n_r <= 0 or (n_t == n_r == 1):
        print("\n No Combining is needed!!")
        return
    n = n_t * n_r
    n_t_tild, n_r_tild = 2 ** math.ceil(math.log2(n_t)), 2 ** math.ceil(math.log2(n_r)),
    n_tild = n_t_tild * n_r_tild
    print(" ================================= Target ==================================")
    print("Nt = " + str(n_t))
    print("Nr = " + str(n_r))
    print("Target QSTBC:")
    _, punctured_indices = produce_reordering_transformation(n_t, n_r)
    target_qstbc = produce_punctured_ea_qstbc(n_tild, n,  punctured_indices)
    print(target_qstbc)
    print("Desired Stream:")
    h_hat = produce_equivalent_channel_vector(n_t * n_r)
    s_hat = np.matmul(target_qstbc, h_hat)
    print(SmartStr(chr(0x015D) + " = ") + SmartStr(target_qstbc) + SmartStr(h_hat) + "  =  " + SmartStr(s_hat))
    print("\n ================================= Result ==================================")
    print("Transmission Matrix:")
    trans_matrix = produce_transmission(n_t, n_r)
    print(trans_matrix)
    print("Combining Scheme:")
    comb_scheme = produce_combining_matrix(n_t, n_r)
    print(comb_scheme)
    print("Channel coefficients transformation:")
    h_transform, _ = produce_reordering_transformation(n_t, n_r)
    print(SmartStr(h_hat) + " = " + SmartStr(h_transform))
    print("\n =============================== Check Result ==============================")
    print("Received Signals Matrix:")
    H = build_channel_matrix(n_t, n_r)
    r = np.matmul(trans_matrix, H)
    print(SmartStr("\nR    =    ") + SmartStr(pandas.DataFrame(r.view(np.ndarray)).to_string()))
    print("Received Stream (after combining):")
    s = do_combining(r)
    print(SmartStr("s = ") + SmartStr(s))
    print("Compare 'Desired Stream' and 'Received Stream':")
    s_hat_to_compare = np.matmul(target_qstbc, h_transform)
    result = (s == s_hat_to_compare)
    print(result)
    if any(result):
        print("\n ========================== Comparison Succeeded !!! ==========================")
    else:
        print("\n ========================== Comparison Failed ===========================")


def present_const_matrices(n_t, n_r):
    if n_t <= 0 or n_r <= 0 or (n_t == n_r == 1):
        print("\n No Combining is needed!!")
        return
    print("Nt = " + str(n_t))
    print("Nr = " + str(n_r))
    print("\n =============================== Transmission Matrices ==============================")
    trans_matrix = produce_transmission(n_t, n_r)
    trans_matrices = transmission_column_representation(trans_matrix)
    print("Each matrix is multiplied from right by the information symbols vector:")
    print(print_information_symbol_vector_str(int(trans_matrices[0].shape[1] / 2)))
    trans_str = SmartStr('')
    for ind, mat in enumerate(trans_matrices):
        trans_str += SmartStr('Tx ' + str(ind) + ':\n' + str(mat) + '    ')
    print(trans_str)
    print("\n =============================== Combining Matrices ==============================")
    comb_matrix = produce_combining_matrix(n_t, n_r)
    comb_matrices = combining_matrix_representation(comb_matrix)
    print("The stream received in rx i (and its conjugates terms):")
    print(print_received_stream_rx_i_str(int(comb_matrices[0].shape[1] / 2)))
    print("The stream received in rx i (and its conjugates terms) is multiplied\n"
          " from left by the corresponding combining matrix and then the results are summed up.")
    comb_str = SmartStr('')
    for ind, mat in enumerate(comb_matrices):
        comb_str += SmartStr('Rx ' + str(ind) + ':\n' + str(mat) + '    ')
    print(comb_str)
