from FunctionsUsingTools import *
from sympy import *
import sys
import codecs

p, q = 1,2


run_title = "4x2"
sys.stdout = codecs.open('Results\\ combining_' + run_title + ".txt", 'w', "utf-8")
print("~~~~~~~~~~~" + run_title + ":: Combining Scheme" + "~~~~~~~~~~~")
present_qstbc_scheme(p, q)
sys.stdout = codecs.open('Results\\ matrices_' + run_title + ".txt", 'w', "utf-8")
print("~~~~~~~~~~~" + run_title + ":: Constant Matrices" + "~~~~~~~~~~~")
present_const_matrices(p, q)
sys.stdout = sys.__stdout__

