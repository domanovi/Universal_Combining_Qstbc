from FunctionsForEA import *
from sympy import *
import sys
import codecs

n_t, n_r = 5, 2


run_title = str(n_t) + "x" + str(n_r)
# sys.stdout = codecs.open('Results\\ combining_' + run_title + ".txt", 'w', "utf-8")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " + run_title + ":: Combining Scheme" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
present_qstbc_scheme(n_t, n_r)
# sys.stdout = codecs.open('Results\\ matrices_' + run_title + ".txt", 'w', "utf-8")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " + run_title + ":: Constant Matrices" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
present_const_matrices(n_t, n_r)
sys.stdout = sys.__stdout__




