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
        return Matrix([["x_" + str(indices[0]), "x_" + str(indices[1])],
                       [conjugate("x_" + str(indices[1])), conjugate("-x_" + str(indices[0]))]])
    else:
        a = produce_abba_qstbc(indices[0: int(len_indices / 2)])
        b = produce_abba_qstbc(indices[int(len_indices / 2):])
        return Matrix([[a, b], [b, a]])


def produce_abba_qstbc_of_order(n):
    """calculates a matrix of n information variables
     forming an ABBA-QSTBC structure"""
    return produce_abba_qstbc(range(n))


def produce_desired_stream_of_order(n):
    """produces the desired stream in ABBA-QSTBC MIMO system with n
    transmitting antennas"""
    c = produce_abba_qstbc_of_order(n)
    h = Matrix(["h_" + str(ind) for ind in range(n)])
    return c * h


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
            result = Matrix([[Symbol("h_" + str(i) + "_-_" + str(j)) for j in rec_indices] for i in range(tx_num)])
            result[:, 1] = conjugate(result[:, 1])
            return result.reshape(tx_num * 2, 1)
        # For more than 2 receivers
        else:
            h_hat_low = produce_equivalent_channel_vec(p, q-1, rec_indices[0: int(rx_num / 2)])
            h_hat_high = produce_equivalent_channel_vec(p, q - 1, rec_indices[int(rx_num / 2):])
            return Matrix([[h_hat_low], [h_hat_high]])


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
            return Matrix([["x_" + str(indices[0]), ],
                           [conjugate("x_" + str(indices[1]))]])
        else:
            c_low = produce_transmission(p-1, q, indices[0: int(len_indices / 2)])
            c_high = produce_transmission(p-1, q, indices[int(len_indices / 2):])
            return Matrix([[c_low, c_high], [c_high, c_low]])
    # For more than 2 receivers
    else:
        c_low = produce_transmission(p, q-1, indices[0: int(len_indices / 2)])
        c_high = produce_transmission(p, q-1, indices[int(len_indices / 2):])
        return Matrix([[c_low], [c_high]])


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
    return Matrix([[Symbol("h_" + str(i) + "_-_" + str(j)) for j in range(rx_num)] for i in range(tx_num)])


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
                result[i, j] = conjugate(mat[i, j])
            elif str(op[i, j]) == "-*":
                result[i, j] = - conjugate(mat[i, j])
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
        result = Matrix([zeros(row_num, 1)])
        for i in range(row_num):
            if i % 2 == 0:
                result[i, 0] = mat[i, 0] + mat[i+1, 1]
            else:
                result[i, 0] = mat[i, 0] + mat[i-1, 1]
        return result
    # For more than 2 receivers
    else:
        a = mat[0: int(row_num/2), 0: int(col_num/2)]
        b = mat[0: int(row_num/2), int(col_num/2):]
        c = mat[int(row_num/2):, 0:int(col_num/2)]
        d = mat[int(row_num/2):, int(col_num/2):]
        return Matrix([[cross_addition(a) + cross_addition(d)],
                       [cross_addition(c) + cross_addition(b)]])


# Printing
# ===================================================================================================
def present_qstbc_scheme(p, q):
    init_printing(use_latex=True)
    print("Operation Matrix:")
    op_mat = produce_operation_matrix(p, q)
    my_print(op_mat)
    print("Transmission Matrix:")
    m = produce_transmission(p, q)
    my_print(m)
    H = build_channel_matrix(p, q)
    s = m*H
    s_hat = operate_matrix(op_mat, s)
    result = cross_addition(s_hat)
    my_print("received stream after combining:")
    pprint(result)
    #
    h_hat = produce_equivalent_channel_vec(p, q)
    stbc = produce_abba_qstbc_of_order(2**p * 2**q)
    y = stbc * h_hat
    print("desired stream:")
    my_print(y)
    print("Desired stream == Received stream?")
    print(y == result)


from sympy import Symbol
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.pretty.stringpict import prettyForm, stringPict
from sympy.printing.pretty.pretty_symbology import _xsym
from sympy.core.function import _coeff_isneg
from sympy.core.function import UndefinedFunction, Function


SUPERSCRIPT_STAR = 0x20F0
BULLET = 0x2022
# override multiplication notation
_xsym['*'] = ('*', ' ' + chr(BULLET) + ' ')


class MyPrettyPrinter(PrettyPrinter):
    def _print_conjugate(self, e):
        """Print conjugate of a symbol with subscript star."""
        sym = e.args
        if not isinstance(*sym, Symbol):
            return super()._print_conjugate(e)
        return prettyForm("{}{}".format(self._print(*sym), chr(SUPERSCRIPT_STAR)))

    def _print_Add(self, expr, order=None):
        if self.order == 'none':
            terms = list(expr.args)
        else:
            terms = self._as_ordered_terms(expr, order=order)
        pforms, indices = [], []

        def pretty_negative(pform, index):
            """Prepend a minus sign to a pretty form. """
            #TODO: Move this code to prettyForm
            if index == 0:
                if pform.height() > 1:
                    pform_neg = '-  '
                else:
                    pform_neg = ' - '
            else:
                pform_neg = '  -  '

            if (pform.binding > prettyForm.NEG
                or pform.binding == prettyForm.ADD):
                p = stringPict(*pform.parens())
            else:
                p = pform
            p = stringPict.next(pform_neg, p)
            # Lower the binding to NEG, even if it was higher. Otherwise, it
            # will print as a + ( - (b)), instead of a - (b).
            return prettyForm(binding=prettyForm.NEG, *p)

        for i, term in enumerate(terms):
            if term.is_Mul and _coeff_isneg(term):
                coeff, other = term.as_coeff_mul(rational=False)
                pform = self._print(Mul(-coeff, *other, evaluate=False))
                pforms.append(pretty_negative(pform, i))
            elif term.is_Rational and term.q > 1:
                pforms.append(None)
                indices.append(i)
            elif term.is_Number and term < 0:
                pform = self._print(-term)
                pforms.append(pretty_negative(pform, i))
            elif term.is_Relational:
                pforms.append(prettyForm(*self._print(term).parens()))
            else:
                pforms.append(self._print(term))

        if indices:
            large = True

            for pform in pforms:
                if pform is not None and pform.height() > 1:
                    break
            else:
                large = False

            for i in indices:
                term, negative = terms[i], False

                if term < 0:
                    term, negative = -term, True

                if large:
                    pform = prettyForm(str(term.p))/prettyForm(str(term.q))
                else:
                    pform = self._print(term)

                if negative:
                    pform = pretty_negative(pform, i)

                pforms[i] = pform

        return MyPrettyForm.__add__(*pforms)


class MyPrettyForm(prettyForm):
    def __add__(self, *others):
        arg = self
        if arg.binding > prettyForm.NEG:
            arg = stringPict(*arg.parens())
        result = [arg]
        for arg in others:
            #add parentheses for weak binders
            if arg.binding > prettyForm.NEG:
                arg = stringPict(*arg.parens())
            #use existing minus sign if available
            if arg.binding != prettyForm.NEG:
                result.append('  +  ')
            result.append(arg)
        return prettyForm(binding=prettyForm.ADD, *stringPict.next(*result))


def my_print(expr):
    print(MyPrettyPrinter().doprint(expr))