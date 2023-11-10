from functions import *

'''
Main file
Newton Rhapson:
    - builds matrices and iterates until convergence
        
'''

def main():
    # Asignment 1: Newton Raphson
    conv_crit = 0.000001
    filenameNR = '/Users/gracedepietro/Desktop/4575algo/PowerFlow_4575/ex_nr.xlsx'
    type1 = 'NR'
    #power_flow_nr_decoup(conv_crit, filenameNR, type1)

    ###################################################################################
    # Assignment 2: Continuation Power Flow
    filenameCPF = '/Users/gracedepietro/Desktop/4575algo/PowerFlow_4575/ex_nr.xlsx'
    contpar = 'volt' #'load'
    beta = np.array([0.5, 0.5])
    alpha = np.array([0.0, 0.0])
    step = 1.0
    #continuousPF(conv_crit, filenameCPF, 'load', alpha, beta, step)

    ###################################################################################
    # Assignment 3: Decoupled Power Flow
    filenameDecoup = '/Users/gracedepietro/Desktop/4575algo/PowerFlow_4575/ex_nr.xlsx'
    # Question 1 included in function "setInitGuess"
    # Question 2
    # 2a
    #primalFDPF(conv_crit, filenameDecoup)
    # 2b
    #dualFastDecoupled(0.0001, filenameDecoup)
    # 2c
    typeSD = 'S_D'
    #power_flow_nr_decoup(conv_crit, filenameDecoup, typeSD)

    ###################################################################################
    # Assignment 4: Distribution Factors and IMML
    # Part 1: Distribution Factors
    filenameDF = '/Users/gracedepietro/Desktop/4575algo/PowerFlow_4575/ex_ass4.xlsx'





if __name__ == "__main__":
    main()