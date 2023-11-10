# Functions file for 4575
# Author: Grace DePietro 2023
from classes import *
import numpy as np
import pandas as pd
import math
import cmath


'''
Function: load file
load initial excel file to get needed lists/info for the rest of the code
'''
def loadFile(filename):
    #read excel file
    #filename = '/Users/gracedepietro/Desktop/4205/project/PowerFlow/' + filename
    initial = pd.read_excel(filename, sheet_name='initial', index_col='bus_num')
    #extract number of buses
    busnum = len(initial.index)
    #extract v and theta values and sub in 0 or 1 for first guess
    v_list = initial.loc[:, 'V'].to_numpy()
    v_list[np.isnan(v_list)] = 1
    v_list = np.array(v_list)
    v_list = v_list.astype('float64')
    t_list = initial.loc[:, 'T'].to_numpy()
    t_list[np.isnan(t_list)] = 0
    t_list = np.array(t_list)
    t_list = t_list.astype('float64')
    #extract p and q list
    #adds NANs to spots where there is no initial
    p_list = initial.loc[:, 'P'].to_numpy()
    q_list = initial.loc[:, 'Q'].to_numpy()


    numP = initial.loc[:, 'P'].count()
    numQ = initial.loc[:, 'Q'].count()
    numT = numP
    numV = numQ
    knownnum = numP + numQ
    line_z = pd.read_excel(filename, sheet_name='line_imp')
    lines = line_z.loc[:, 'line'].to_numpy()
    r_list = line_z.loc[:, 'R'].to_numpy()
    r_list = r_list.astype('float64')
    x_list = line_z.loc[:, 'X'].to_numpy()
    x_list = x_list.astype('float64')

    return v_list, t_list, p_list, q_list, lines, r_list, x_list, knownnum, busnum, line_z, numT, numV

"""
Function: get initial matrices
Gets values from the lists from excel and puts them in the mismatch vector and unknown vector
Parameters:
    busnum - total number of buses
    xmat - empty matrix of size busnum rows and 1 column (unknowns matrix)
    knowns - empty matrix of size busnum rows and 1 column (knowns matrix)
    p_list - p vals
    q_list - q vals
Returns:
    knowns - the filled known matrix with name and values
    xmat - the filled unknown matrix with name only

"""
def getInitMats(xmat, knowns, p_list, q_list, busnum):
    sumcount = 0
    for i in range(busnum):
        # check if there is a value for initial P
        # add a T to the xmat and a P to the knowns
        if np.isnan(p_list[i]) == False:
            xmat[sumcount].name = "T" + str(i + 1)
            knowns[sumcount].name = "P" + str(i + 1)
            knowns[sumcount].val = p_list[i]
            sumcount += 1
    newcount = 0
    for j in range(busnum):
        # check if there is a value for initial Q
        # add a V to the xmatrix and a Q to the knowns
        if np.isnan(q_list[j]) == False:
            xmat[newcount + sumcount].name = "V" + str(j + 1)
            knowns[newcount + sumcount].name = "Q" + str(j + 1)
            knowns[newcount + sumcount].val = q_list[j]
            newcount += 1


"""
Function: Set initial guess
This function sets the initial guesses of 1.0pu for voltage and 0 for angle
Parameters:
    knownnum - number of knowns in the system
    xmat - matrix of unknowns with only names
Returns
    xmat - matrix of unknowns with names and initial guesses
"""
def setInitGuess(knownNum, xmat, v_list, t_list):
    for i in range(int(knownNum)):
        if "V" in xmat[i].name:
            xmat[i].val = v_list[int(xmat[i].name[1])-1]
        elif "T" in xmat[i].name:
            xmat[i].val = t_list[int(xmat[i].name[1])-1]


"""
Function: print matrix
This will print a known or unknown matrix (rows amount = busnum, column = 1)
Parameters:
    num - number of elem for rows in matrix
    xmat - matrix
Returns:
    nothing, just print
"""
def printMat(num, xmat):
    for i in range(int(num)):
        print(xmat[i].name + ": " + str(xmat[i].val))


"""
Function: print multi matrix
This will print a matrix (probably Ybus, jacobian)
Parameters:
    num - number of elem for rows and columns in matrix
    xmat - matrix
Returns:
    nothing, just print
"""
def printMultiMat(num, mat, jac):
    str_print = ['' for i in range(int(num))]
    for i in range(int(num)):
        for j in range(int(num)):
            if jac == False:
                #print(mat[i][j].name + ", " + str(mat[i][j].val))
                str_print[i] += mat[i][j].name + ": " + str(mat[i][j].val)
            else:
                #print(mat[i][j].name + ", " + str(mat[i][j].val) + ", " + str(mat[i][j].type))
                str_print[i] += mat[i][j].name + ": " + str(mat[i][j].val)
            if j != int(num):
                str_print[i] += ', '
        print(str_print[i])


"""
Function: name Y bus
This will set the names of the Ybus matrix things.

Parameters:
    busnum - number of nuses for matrix size indeces
    yBus - matrix with yBus names and vals
"""

def nameYbus(busnum, yBus):
    for i in range(int(busnum)):
        for j in range(int(busnum)):
            yBus[i][j].name = "Y" + str(i + 1) + str(j + 1)



'''
Function: pi line model
for creating a 2x2 admittance matrix for a one line system. Part 1 of the Cutsem method.
'''
def piLine(r_list, x_list, y_mini, lines):
    line_num = lines.size
    #shape: 0 is 11, 1 is 12, 2 is 21, and 3 is 22
    for i in range(line_num):
        y_mini[i][0] = 1 / complex(r_list[i], x_list[i])
        y_mini[i][1] = -1/complex(r_list[i],x_list[i])
        y_mini[i][2] = y_mini[i][1]
        y_mini[i][3] = y_mini[i][0]


'''
Function: Ybus calculations with cutsem algorithm/pi line model
'''
def yBusCutsemCalc(y_mini, lines, yBus):
    for i in range(len(lines)):
        str_temp = str(lines[i])
        a = int(str_temp[0])-1
        b = int(str_temp[1])-1
        yBus[a][a].val += y_mini[i][0]
        yBus[a][b].val += y_mini[i][1]
        yBus[b][a].val += y_mini[i][2]
        yBus[b][b].val += y_mini[i][3]
'''
Function: calculate uij
'''
def uij(gij, bij, thetai, thetaj):
    return (-gij * np.sin(thetai - thetaj) - (-bij * np.cos(thetai - thetaj)))


'''
Function: calculate tij
'''
def tij(gij, bij, thetai, thetaj):
    return (-gij * np.cos(thetai - thetaj) + (-bij * np.sin(thetai - thetaj)))


'''
Function: calculate P value
'''
def calcP(i, yBus, busnum, T, V):
    p = 0
    for j in range(busnum):
        p += V[i]*V[j]*abs(yBus[i][j].val)*np.cos(T[i]-T[j]-cmath.phase(yBus[i][j].val))
    return p

'''
Function: calculate Q value
'''
def calcQ(i, yBus, busnum, T, V):
    q = 0
    for j in range(busnum):
        q += V[i]*V[j]*abs(yBus[i][j].val)*np.sin(T[i]-T[j]-cmath.phase(yBus[i][j].val))
    return q


'''
Function: calculate partial derivative dPi / dTi
'''
def dpidti(i, V, yBus, T, busnum):
    i = int(i)
    sum = 0
    for j in range(busnum):
        if j != i:
            sum += V[j] * uij(yBus[i][j].val.real, yBus[i][j].val.imag, T[i], T[j])
    return sum * V[i]


'''
Function: calculate partial derivative dPi / dTj
'''
def dpidtj(i, j, V, yBus, T):
    i = int(i)
    j = int(j)
    return -V[i] * V[j] * uij(yBus[i][j].val.real, yBus[i][j].val.imag, T[i], T[j])


'''
Function: calculate partial derivative dPi / dVi
'''
def dpidvi(i, V, yBus, T, busnum):
    i = int(i)
    sum = 2 * V[i] * yBus[i][i].val.real
    for j in range(busnum):
        if j != i:
            sum += -V[j] * tij(yBus[i][j].val.real, yBus[i][j].val.imag, T[i], T[j])
    return sum


'''
Function: calculate partial derivative dPi / dQj
'''
def dpidvj(i, j, V, yBus, T):
    i = int(i)
    j = int(j)
    return -V[i] * tij(yBus[i][j].val.real, yBus[i][j].val.imag, T[i], T[j])


'''
Function: calculate partial derivative dQi / dTi

'''
def dqidti(i, V, yBus, T, busnum):
    i = int(i)
    sum = 0
    for j in range(busnum):
        if j != i:
            sum += V[j] * tij(yBus[i][j].val.real, yBus[i][j].val.imag, T[i], T[j])
    return sum * -V[i]


'''
Function: calculate partial derivative dQi / dTj

'''
def dqidtj(i, j, V, yBus, T):
    i = int(i)
    j = int(j)
    return V[i] * V[j] * tij(yBus[i][j].val.real, yBus[i][j].val.imag, T[i], T[j])


'''
Function: calculate partial derivative dQi / dVi
'''
def dqidvi(i, V, yBus, T, busnum):
    i = int(i)
    sum = -2 * V[i] * yBus[i][i].val.imag
    for j in range(busnum):
        if j != i:
            sum += -V[j] * uij(yBus[i][j].val.real, yBus[i][j].val.imag, T[i], T[j])
    return sum


'''
Function: calculate partial derivative dQi / dVj
'''
def dqidvj(i, j, V, yBus, T):
    i = int(i)
    j = int(j)
    x = (-V[i] * uij(yBus[i][j].val.real, yBus[i][j].val.imag, T[i], T[j]))
    return x

'''
Function: calculate partial derivative dPi / dTi for dual FDLF
'''
def dpidtiDFDLF(i, v_list, yBus):
    return -v_list[i]*yBus[i][i].val.imag

'''
Function: calculate partial derivative dPi / dTj for dual FDLF
'''
def dpidtjDFDLF(i, j, v_list, yBus):
    return -v_list[i]*yBus[i][j].val.imag

'''
Function: calculate partial derivative dQi / dVi for primal FDLF
'''
def dqidviPFDLF(i, v_list, yBus):
    return -v_list[i]*yBus[i][i].val.imag

'''
Function: calculate partial derivative dQi / dVj for primal FDLF
'''
def dqidvjPFDLF(i, j, v_list, yBus):
    return -v_list[i]*yBus[i][j].val.imag

'''
Function: name jacobian element
names the elements of the jacobian matrix for easier readability/understanding
'''
def nameJacElem(knownnum, knowns, xmat, jacobian):
    for i in range(knownnum):
        for j in range(knownnum):
            # Ps
            if knowns[i].name[0] == 'P' and xmat[j].name[0] == 'T':
                jacobian[i][j].name = 'dp' + str(knowns[i].name[1]) + 'dt' + str(xmat[j].name[1])
                if knowns[i].name[1] == xmat[j].name[1]:
                    jacobian[i][j].type = 'dpidti'
                else:
                    jacobian[i][j].type = 'dpidtj'
            if knowns[i].name[0] == 'P' and xmat[j].name[0] == 'V':
                jacobian[i][j].name = 'dp' + str(knowns[i].name[1]) + 'dv' + str(xmat[j].name[1])
                if knowns[i].name[1] == xmat[j].name[1]:
                    jacobian[i][j].type = 'dpidvi'
                else:
                    jacobian[i][j].type = 'dpidvj'
            # Qs
            if knowns[i].name[0] == 'Q' and xmat[j].name[0] == 'T':
                jacobian[i][j].name = 'dq' + str(knowns[i].name[1]) + 'dt' + str(xmat[j].name[1])
                if knowns[i].name[1] == xmat[j].name[1]:
                    jacobian[i][j].type = 'dqidti'
                else:
                    jacobian[i][j].type = 'dqidtj'
            if knowns[i].name[0] == 'Q' and xmat[j].name[0] == 'V':
                jacobian[i][j].name = 'dq' + str(knowns[i].name[1]) + 'dv' + str(xmat[j].name[1])
                if knowns[i].name[1] == xmat[j].name[1]:
                    jacobian[i][j].type = 'dqidvi'
                else:
                    jacobian[i][j].type = 'dqidvj'


'''
Function: calculate jacobian elements and update matrix
'''
def calcJacElems(knownnum, jacobian, ybus, t_list, v_list, busnum):
    for i in range(knownnum):
        for j in range(knownnum):
            # Ps
            # this i and j isnt from the loop. its from the value from P/Q and V/T
            i_temp = int(jacobian[i][j].name[2])-1
            j_temp = int(jacobian[i][j].name[5])-1
            if jacobian[i][j].type == 'dpidti':
                jacobian[i][j].val = dpidti(i_temp, v_list, ybus, t_list, busnum)
            elif jacobian[i][j].type == 'dpidtj':
                jacobian[i][j].val = dpidtj(i_temp, j_temp, v_list, ybus, t_list)
            elif jacobian[i][j].type == 'dpidvi':
                jacobian[i][j].val = dpidvi(i_temp, v_list, ybus, t_list, busnum)
            elif jacobian[i][j].type == 'dpidvj':
                jacobian[i][j].val = dpidvj(i_temp, j_temp, v_list, ybus, t_list)
            elif jacobian[i][j].type == 'dqidti':
                jacobian[i][j].val = dqidti(i_temp, v_list, ybus, t_list, busnum)
            elif jacobian[i][j].type == 'dqidtj':
                jacobian[i][j].val = dqidtj(i_temp, j_temp, v_list, ybus, t_list)
            elif jacobian[i][j].type == 'dqidvi':
                jacobian[i][j].val = dqidvi(i_temp, v_list, ybus, t_list, busnum)
            elif jacobian[i][j].type == 'dqidvj':
                jacobian[i][j].val = dqidvj(i_temp, j_temp, v_list, ybus, t_list)
            else:
                print('error')

'''
Function: calculate jacobian elements and update matrix for Decoupled load flow
sets specific elements to zero
'''
def calcJacElemsDLF(knownnum, jacobian, ybus, t_list, v_list, busnum):
    for i in range(knownnum):
        for j in range(knownnum):
            # Ps
            # this i and j isnt from the loop. its from the value from P/Q and V/T
            i_temp = int(jacobian[i][j].name[2])-1
            j_temp = int(jacobian[i][j].name[5])-1
            if jacobian[i][j].type == 'dpidti':
                jacobian[i][j].val = dpidti(i_temp, v_list, ybus, t_list, busnum)
            elif jacobian[i][j].type == 'dpidtj':
                jacobian[i][j].val = dpidtj(i_temp, j_temp, v_list, ybus, t_list)
            elif jacobian[i][j].type == 'dpidvi':
                jacobian[i][j].val = 0
            elif jacobian[i][j].type == 'dpidvj':
                jacobian[i][j].val = 0
            elif jacobian[i][j].type == 'dqidti':
                jacobian[i][j].val = 0
            elif jacobian[i][j].type == 'dqidtj':
                jacobian[i][j].val = 0
            elif jacobian[i][j].type == 'dqidvi':
                jacobian[i][j].val = dqidvi(i_temp, v_list, ybus, t_list, busnum)
            elif jacobian[i][j].type == 'dqidvj':
                jacobian[i][j].val = dqidvj(i_temp, j_temp, v_list, ybus, t_list)
            else:
                print('error')

'''
Function: Function: calculate jacobian elements and update matrix for Fast Decoupled load flow
sets specific elements to zero
'''

'''
Function: iterate
should update P, Q, known matrix, unknown matrix, and jacobian
Parameters:
    - knownnum: int of # elems in mismatch P/Q vector
    - jacobian: jacobian vector (empty with names only)
    - ybus: admittance matrix
    - t_list: list of thetas
    - v_lists: list of voltage magnitudes
    - knowns: mismatch P/Q vector
    - xmat: unknowns V/T vector
    - busnum: int total # buses
    - type for "NR" (newton rhapson) or 'S_D' (standard decoupled)
'''
def iterate(knownnum, jacobian, ybus, t_list, v_list, knowns, xmat, busnum, type):
    #first calculate the jacobian matrix
    if type == 'NR':
        calcJacElems(knownnum, jacobian, ybus, t_list, v_list, busnum)
    elif type == 'S_D':
        calcJacElemsDLF(knownnum, jacobian, ybus, t_list, v_list, busnum)
    else:
        print('error thrown in which type of jacobian calculation')

    #make temp knowns matrix without the names
    new_knowns = [0 for i in range(knownnum)]
    net_injections = [0 for i in range(knownnum)]
    for i in range(knownnum):
        new_knowns[i] = knowns[i].val
    for i in range(knownnum):
        #for each known value, calculate the new value of P or Q and subtract them
        num = int(knowns[i].name[1])-1
        type = knowns[i].name[0]
        if type == 'P':
            #Note: change generating/not +/- for P and Q IN EXCEL
            new_p = calcP(num, ybus, busnum, t_list, v_list)
            net_injections[i] = new_p
            new_knowns[i] = new_knowns[i] - new_p
        else:
            #Note: change generating/not +/- for P and Q IN EXCEL
            new_q = calcQ(num, ybus, busnum, t_list, v_list)
            net_injections[i] = new_q
            new_knowns[i] = new_knowns[i] - new_q
    print("Net Injections: ")
    for i in range(knownnum):
        print(knowns[i].name, ': ', net_injections[i])


    temp_jac = [[0 for i in range(int(knownnum))] for j in range(int(knownnum))]
    for i in range(knownnum):
        for j in range(knownnum):
            temp_jac[i][j] = jacobian[i][j].val
    corrections = np.linalg.solve(temp_jac, new_knowns)
    for j in range(knownnum):
        xmat[j].val += corrections[j]
        temp_num = int(xmat[j].name[1])
        temp_val = xmat[j].val
        if xmat[j].name[0] == "T":
            t_list[(temp_num-1)] = temp_val
        elif xmat[j].name[0] == "V":
            v_list[(temp_num-1)] = temp_val
        else:
            print("error thrown in updating v and t lists")


    print("Jacobian: ")
    printMultiMat(knownnum, jacobian, True)
    print("Corrections Vector (dV, dT): ")
    print(corrections)
    print("RHS Vector (dP, dQ): ")
    print(new_knowns)
    print("New voltage magnitudes and angles (in degrees):")
    for i in range(busnum):
        print("|V|", i + 1, ": ", "{:.4f}".format(v_list[i]), "\t\t", "Theta", i + 1, ": ",
              "{:.4f}".format(t_list[i] * 180 / math.pi))

    return corrections

"""
Function: loop normal
loop until convergence for normal newton raphson
parameters:
    - knowns: known matrix mismatch vector
    - knownnum: length of known matrix
    - jacobian: jacobian matrix
    - yBus: admittance matrix
    - t_list: list of angles
    - v_list: list of voltage mags
    - xmat: unkown vector of Vs and Ts
    - busnum: total number of buses
    - conv_crit: convergence criteria (small number)
    - type: string of "NR" (newton raphson) or "S_D" (standard decoupled)
"""
def loop_normal(knowns, knownnum, jacobian, yBus, t_list, v_list, xmat, busnum, conv_crit, type):
    convergence = False
    itno = 0
    while not (convergence or itno > 10):
        itno += 1
        print("\n\nIteration #" + str(itno))
        corrections = iterate(knownnum, jacobian, yBus, t_list, v_list, knowns, xmat, busnum, type)

        # Check corrections matrix for convergence
        count = 0
        for i in range(corrections.size):
            if abs(corrections[i]) > conv_crit:
                cur = abs(corrections[i])
                count += 1

        convergence = count == 0



"""
Function: power flow in newton raphson and decoupled
This function runs either newton raphson or decoupled power flow based on the type parameter entered
Parameters:
    - conv_crit: float of the convergence criteria
    - filename: string excel file path name with input data
    - type: string of "NR" (newton raphson) or "S_D" (standard decoupled)
"""
def power_flow_nr_decoup(conv_crit, filename, type):
    # load info from file, set up matrices:
    stuff = loadFile(filename)
    v_list = stuff[0]
    t_list = stuff[1]
    p_list = stuff[2]
    q_list = stuff[3]
    lines = stuff[4]
    r_list = stuff[5]
    x_list = stuff[6]
    knownnum = stuff[7]
    busnum = stuff[8]
    line_z = stuff[9]
    numT = stuff[10]
    numV = stuff[11]

    knowns = [VarMat() for i in range(int(knownnum))]
    xmat = [VarMat() for j in range(int(knownnum))]


    getInitMats(xmat, knowns, p_list, q_list, busnum)
    setInitGuess(knownnum, xmat, v_list, t_list)


    yBus = [[VarMat() for i in range(int(busnum))] for j in range(int(busnum))]

    nameYbus(busnum, yBus)


    y_mini = [[complex(0, 0) for i in range(4)] for j in range(lines.size)]
    piLine(r_list, x_list, y_mini, lines)
    yBusCutsemCalc(y_mini, lines, yBus)
    #printMultiMat(busnum, yBus, False)
    jacobian = [[JacElem() for i in range(int(knownnum))] for j in range(int(knownnum))]
    nameJacElem(knownnum, knowns, xmat, jacobian)

    # run loop until convergence criteria is met
    # jacobian type will change for newton raphson or decoupled
    loop_normal(knowns, knownnum, jacobian, yBus, t_list, v_list, xmat, busnum, conv_crit, type)


    for i in range(busnum):
        if np.isnan(p_list[i]):
            p_list[i] = calcP(i, yBus, busnum, t_list, v_list) # calcPVal(i, yBus, busnum, t_list, v_list)
        if np.isnan(q_list[i]):
            q_list[i] = calcQ(i, yBus, busnum, t_list, v_list) # calcQVal(i, yBus, busnum, t_list, v_list)
    print('Final P and Q Values: ')
    for i in range(busnum):
        print("P", i + 1, ": ", "{:.4f}".format(p_list[i]), "\t\t\t", "Q", i + 1, ": ", "{:.4f}".format(q_list[i]))
    print('Final V and T Values: ')
    for i in range(busnum):
        print("V", i + 1, ": ", "{:.4f}".format(v_list[i]), "\t\t\t", "T", i + 1, ": ", "{:.4f}".format(t_list[i]))

"""
Function: extend jacobian for continuous power flow
takes in the normal and extends it. adds a zero in the last row based on the parameter "where" (int)(column no zero based)
"""
def extendJacCPF(knownnum, jacobian_new, jacobian, alpha, beta, where):
    numhalf = int(knownnum/2)
    for i in range(knownnum+1):
        for j in range(knownnum+1):
            if i < knownnum and j < knownnum: # for 3 bus: rows 1-4 and columns 1-4
                jacobian_new[i][j] = jacobian[i][j].val
            elif i < knownnum: # for 3 bus: rows 1-4
                if j == knownnum: # last column (5)
                    if i < numhalf: # if first or second row
                        jacobian_new[i][j] = beta[i]
                    else: # if third or 4th row
                        intcur = i-numhalf
                        jacobian_new[i][j] = alpha[intcur]
            else: # fifth row
                if j != where: # columns 1-4
                    jacobian_new[i][j] = 0.0
                else: # last column
                    jacobian_new[i][j] = 1.0


def continuousPF(conv_crit, filename, contpar, alpha, beta, step):
    # load info from file, set up matrices:
    stuff = loadFile(filename)
    v_list = stuff[0]
    t_list = stuff[1]
    p_list = stuff[2]
    q_list = stuff[3]
    lines = stuff[4]
    r_list = stuff[5]
    x_list = stuff[6]
    knownnum = stuff[7]
    busnum = stuff[8]
    line_z = stuff[9]
    numT = stuff[10]
    numV = stuff[11]

    knowns = [VarMat() for i in range(int(knownnum))]
    xmat = [VarMat() for j in range(int(knownnum))]

    getInitMats(xmat, knowns, p_list, q_list, busnum)
    setInitGuess(knownnum, xmat, v_list, t_list)

    yBus = [[VarMat() for i in range(int(busnum))] for j in range(int(busnum))]

    nameYbus(busnum, yBus)

    y_mini = [[complex(0, 0) for i in range(4)] for j in range(lines.size)]
    piLine(r_list, x_list, y_mini, lines)
    yBusCutsemCalc(y_mini, lines, yBus)
    # printMultiMat(busnum, yBus, False)
    jacobian = [[JacElem() for i in range(int(knownnum))] for j in range(int(knownnum))]
    nameJacElem(knownnum, knowns, xmat, jacobian)

    # Perform one Newton-Raphson iteration
    # updates t_list and v_list
    nr = 0
    while nr < 4:
        corrections = iterate(knownnum, jacobian, yBus, t_list, v_list, knowns, xmat, busnum, 'NR')
        nr += 1

    # Perform one iteration with a modified jacobian

    jacobian_new = [[0 for i in range(int(knownnum+1))] for j in range(int(knownnum+1))]
    extendJacCPF(knownnum, jacobian_new, jacobian, alpha, beta, knownnum-1)

    b = np.array([0.0, 0.0, 0.0, 0.0, 1.0])
    # solve the matrices
    corrections = np.linalg.solve(jacobian_new, b)

    # update state variables using step size
    print("Sensitivities: ", corrections)
    for i in range(knownnum):
        if i < numT:
            t_list[i] += step*corrections[i]
            p_list[i] += step*beta[i]
        else:
            v_list[i-numT] += step*corrections[i]
            q_list[i-numT] += step*alpha[i-numT]
    print("t list: ", t_list)
    print("v list: ", v_list)
    print("p list: ", p_list)
    print("q list: ", q_list)
    p01 = p_list[0]
    p02 = p_list[1]
    q01 = q_list[0]
    q02 = q_list[1]

    if contpar == 'load':
        itno = 0
        convergence = False
        while not (convergence or itno > 10):
            # calculate elements of jacobia
            calcJacElems(knownnum, jacobian, yBus, t_list, v_list, busnum)
            # extend jacobian with alpha and beta
            jacobian_new = [[0 for i in range(int(knownnum + 1))] for j in range(int(knownnum + 1))]
            extendJacCPF(knownnum, jacobian_new, jacobian, alpha, beta, knownnum-1)
            # calc net injections
            p1 = calcP(0, yBus, busnum, t_list, v_list)
            p2 = calcP(1, yBus, busnum, t_list, v_list)
            q1 = calcQ(0, yBus, busnum, t_list, v_list)
            q2 = calcQ(1, yBus, busnum, t_list, v_list)
            # create dp dq vector
            b = np.array([p01-p1, p02-p2, q01-q1, q02-q2, 0])
            corrections = np.linalg.solve(jacobian_new, b)
            # update v and t lists
            for i in range(knownnum):
                if i < numT:
                    t_list[i] += corrections[i]
                else:
                    v_list[i - numT] += corrections[i]

            print('Iteration number ', itno, ':')
            print('Net Injections: ')
            print("P", 1, ": ", "{:.4f}".format(p1), "\t\t\t", "Q", 1, ": ", "{:.4f}".format(q1))
            print("P", 2, ": ", "{:.4f}".format(p2), "\t\t\t", "Q", 2, ": ", "{:.4f}".format(q2))

            print("Corrections: ")
            print(corrections)

            print("dP dQ: ")
            print("P", 1, ": ", "{:.4f}".format(p01-p1), "\t\t\t", "Q", 1, ": ", "{:.4f}".format(q01-q1))
            print("P", 2, ": ", "{:.4f}".format(p01-p2), "\t\t\t", "Q", 2, ": ", "{:.4f}".format(q02-q2))
            print('v1: ', v_list[0], '  v2: ', v_list[1], 'theta1: ', t_list[0], '  theta2: ', t_list[1])

            #check for convergence
            count = 0
            for i in range(corrections.size):
                if abs(corrections[i]) > conv_crit:
                    cur = abs(corrections[i])
                    count += 1
            convergence = count == 0
            itno += 1

    elif contpar == 'volt':
        itno = 0
        convergence = False
        while not (convergence or itno > 10):
            # calculate elements of jacobia
            calcJacElems(knownnum, jacobian, yBus, t_list, v_list, busnum)
            # extend jacobian with alpha and beta
            jacobian_new = [[0 for i in range(int(knownnum + 1))] for j in range(int(knownnum + 1))]
            extendJacCPF(knownnum, jacobian_new, jacobian, alpha, beta, knownnum-2)
            # calc net injections
            p1 = calcP(0, yBus, busnum, t_list, v_list)
            p2 = calcP(1, yBus, busnum, t_list, v_list)
            q1 = calcQ(0, yBus, busnum, t_list, v_list)
            q2 = calcQ(1, yBus, busnum, t_list, v_list)
            # create dp dq vector
            b = np.array([p01-p1, p02-p2, q01-q1, q02-q2, 0])
            corrections = np.linalg.solve(jacobian_new, b)
            # update v and t lists
            for i in range(knownnum):
                if i < numT:
                    t_list[i] += corrections[i]
                else:
                    v_list[i - numT] += corrections[i]

            print('Iteration number ', itno, ':')
            print('Net Injections: ')
            print("P", 1, ": ", "{:.4f}".format(p1), "\t\t\t", "Q", 1, ": ", "{:.4f}".format(q1))
            print("P", 2, ": ", "{:.4f}".format(p2), "\t\t\t", "Q", 2, ": ", "{:.4f}".format(q2))

            print("Corrections: ")
            print(corrections)

            print("dP dQ: ")
            print("P", 1, ": ", "{:.4f}".format(p01-p1), "\t\t\t", "Q", 1, ": ", "{:.4f}".format(q01-q1))
            print("P", 2, ": ", "{:.4f}".format(p01-p2), "\t\t\t", "Q", 2, ": ", "{:.4f}".format(q02-q2))
            print('v1: ', v_list[0], '  v2: ', v_list[1], 'theta1: ', t_list[0], '  theta2: ', t_list[1])

            s = corrections[knownnum]
            for i in range(knownnum):
                if i < numT:
                    p_list[i] += s * beta[i]
                else:
                    q_list[i - numT] += s * alpha[i - numT]

            #check for convergence
            count = 0
            for i in range(corrections.size):
                if abs(corrections[i]) > conv_crit:
                    cur = abs(corrections[i])
                    count += 1
            convergence = count == 0
            itno += 1

    else:
        print("entered wrong type of input for continuous power flow calculation, please change string to 'load' or 'volt'")

    #print final values
    for i in range(busnum):
        if np.isnan(p_list[i]):
            p_list[i] = calcP(i, yBus, busnum, t_list, v_list)
        if np.isnan(q_list[i]):
            q_list[i] = calcQ(i, yBus, busnum, t_list, v_list)
    print('Final P and Q Values: ')
    for i in range(busnum):
        print("P", i + 1, ": ", "{:.4f}".format(p_list[i]), "\t\t\t", "Q", i + 1, ": ", "{:.4f}".format(q_list[i]))
    print('Final V and T Values: ')
    for i in range(busnum):
        print("V", i + 1, ": ", "{:.4f}".format(v_list[i]), "\t\t\t", "T", i + 1, ": ", "{:.4f}".format(t_list[i]))



"""
Function: Primal fast decoupled load flow
for assignment 3 question 2a
parameters:
    - conv_crit: float of convergence criteria
    - filename: string of path to excel file
"""
def primalFDPF(conv_crit, filename):
    stuff = loadFile(filename)
    v_list = stuff[0]
    t_list = stuff[1]
    p_list = stuff[2]
    q_list = stuff[3]
    lines = stuff[4]
    r_list = stuff[5]
    x_list = stuff[6]
    knownnum = stuff[7]
    busnum = stuff[8]
    line_z = stuff[9]
    numT = stuff[10]
    numV = stuff[11]

    knowns = [VarMat() for i in range(int(knownnum))]
    xmat = [VarMat() for j in range(int(knownnum))]


    getInitMats(xmat, knowns, p_list, q_list, busnum)
    setInitGuess(knownnum, xmat, v_list, t_list)


    yBus = [[VarMat() for i in range(int(busnum))] for j in range(int(busnum))]

    nameYbus(busnum, yBus)


    y_mini = [[complex(0, 0) for i in range(4)] for j in range(lines.size)]
    piLine(r_list, x_list, y_mini, lines)
    yBusCutsemCalc(y_mini, lines, yBus)
    itno = 0
    convergence = False
    p01 = p_list[0]
    p02 = p_list[1]
    q01 = q_list[0]
    q02 = q_list[1]

    while not (convergence or itno > 10):
        # create "jacobian"
        b_primal = np.matrix([[yBus[0][0].val.imag, yBus[0][1].val.imag], [yBus[0][1].val.imag, yBus[1][1].val.imag]])

        # update P values and build mismatch
        p1 = calcP(0, yBus, busnum, t_list, v_list)
        p2 = calcP(1, yBus, busnum, t_list, v_list)
        b_p = np.array([p01-p1, p02-p2])

        # build change P / V
        dp_v = b_p/[v_list[0], v_list[1]]

        # solve for correction vector of theta
        corrections_p = np.linalg.solve(b_primal, dp_v)

        # update t_list
        t_list[0] -= corrections_p[0]
        t_list[1] -= corrections_p[1]

        # update Q values and build mismatch
        q1 = calcQ(0, yBus, busnum, t_list, v_list)
        q2 = calcQ(1, yBus, busnum, t_list, v_list)
        b_q = np.array([q01-q1, q02-q1])

        # build change Q / V
        dq_v = b_q/[v_list[0], v_list[1]]
        # get corrections vector
        corrections_q = np.linalg.solve(b_primal, dq_v)

        # update v_list
        v_list[0] -= corrections_q[0]
        v_list[1] -= corrections_q[1]

        # combine corrections vector
        corrections = np.array([corrections_p[0], corrections_p[1], corrections_q[0], corrections_q[1]])
        print('Iteration number ', itno, ':')
        print('Net Injections: ')
        print("P", 1, ": ", "{:.4f}".format(p1), "\t\t\t", "Q", 1, ": ", "{:.4f}".format(q1))
        print("P", 2, ": ", "{:.4f}".format(p2), "\t\t\t", "Q", 2, ": ", "{:.4f}".format(q2))

        print("Corrections: ")
        print(corrections)

        print("dP dQ: ")
        print("P", 1, ": ", "{:.4f}".format(p01-p1), "\t\t\t", "Q", 1, ": ", "{:.4f}".format(q01-q1))
        print("P", 2, ": ", "{:.4f}".format(p01-p2), "\t\t\t", "Q", 2, ": ", "{:.4f}".format(q02-q2))

        count = 0
        for i in range(corrections.size):
            if abs(corrections[i]) > conv_crit:
                cur = abs(corrections[i])
                count += 1
        convergence = count == 0
        itno += 1

    for i in range(busnum):
        if np.isnan(p_list[i]):
            p_list[i] = calcP(i, yBus, busnum, t_list, v_list)
        if np.isnan(q_list[i]):
            q_list[i] = calcQ(i, yBus, busnum, t_list, v_list)
    print('Final P and Q Values: ')
    for i in range(busnum):
        print("P", i + 1, ": ", "{:.4f}".format(p_list[i]), "\t\t\t", "Q", i + 1, ": ", "{:.4f}".format(q_list[i]))
    print('Final V and T Values: ')
    for i in range(busnum):
        print("V", i + 1, ": ", "{:.4f}".format(v_list[i]), "\t\t\t", "T", i + 1, ": ", "{:.4f}".format(t_list[i]))


"""
Function: Dual fast decoupled load flow
for assignment 3 question 2b
parameters:
    - conv_crit: float of convergence criteria
    - filename: string of path to excel file
"""
def dualFastDecoupled(conv_crit, filename):

    stuff = loadFile(filename)
    v_list = stuff[0]
    t_list = stuff[1]
    p_list = stuff[2]
    q_list = stuff[3]
    lines = stuff[4]
    r_list = stuff[5]
    x_list = stuff[6]
    knownnum = stuff[7]
    busnum = stuff[8]
    line_z = stuff[9]
    numT = stuff[10]
    numV = stuff[11]

    knowns = [VarMat() for i in range(int(knownnum))]
    xmat = [VarMat() for j in range(int(knownnum))]


    getInitMats(xmat, knowns, p_list, q_list, busnum)
    setInitGuess(knownnum, xmat, v_list, t_list)


    yBus = [[VarMat() for i in range(int(busnum))] for j in range(int(busnum))]

    nameYbus(busnum, yBus)


    y_mini = [[complex(0, 0) for i in range(4)] for j in range(lines.size)]
    piLine(r_list, x_list, y_mini, lines)
    yBusCutsemCalc(y_mini, lines, yBus)
    itno = 0
    convergence = False
    p01 = p_list[0]
    p02 = p_list[1]
    q01 = q_list[0]
    q02 = q_list[1]

    while not (convergence or itno > 10):
        # create "jacobian"
        b_primal = np.matrix([[yBus[0][0].val.imag, yBus[0][1].val.imag], [yBus[0][1].val.imag, yBus[1][1].val.imag]])

        q1 = calcQ(0, yBus, busnum, t_list, v_list)
        q2 = calcQ(1, yBus, busnum, t_list, v_list)

        # Build mismatch for Q
        b_q = np.array([q01-q1 , q02-q2 ])

        # Build deltaQ/V
        deltaQ_V = b_q / [v_list[0], v_list[1]]

        # Calculate correction of v
        corrections_q = np.linalg.solve(b_primal, deltaQ_V)

        # Update v
        v_list[0] -= corrections_q[0]
        v_list[1] -= corrections_q[1]

        # Find P mismatch
        # Update P values
        p1 = calcP(0, yBus, busnum, t_list, v_list)
        p2 = calcP(1, yBus, busnum, t_list, v_list)

        # Build mismatch for P
        b_p = np.array([p01-p1, p02-p2])

        # Build deltaP/V
        deltaP_V = b_p / [v_list[0], v_list[1]]

        # Calculate correction of theta
        corrections_p = np.linalg.solve(b_primal, deltaP_V)

        # Update theta
        t_list[0] -= corrections_p[0]
        t_list[1] -= corrections_p[1]

        # combine corrections vector
        corrections = np.array([corrections_p[0], corrections_p[1], corrections_q[0], corrections_q[1]])
        print('Iteration number', itno+1, ':')
        print('Net Injections: ')
        print("P", 1, ": ", "{:.4f}".format(p1), "\t\t\t", "Q", 1, ": ", "{:.4f}".format(q1))
        print("P", 2, ": ", "{:.4f}".format(p2), "\t\t\t", "Q", 2, ": ", "{:.4f}".format(q2))

        print("Corrections: ")
        print(corrections)

        print("dP dQ: ")
        print("P", 1, ": ", "{:.4f}".format(p01 - p1), "\t\t\t", "Q", 1, ": ", "{:.4f}".format(q01 - q1))
        print("P", 2, ": ", "{:.4f}".format(p01 - p2), "\t\t\t", "Q", 2, ": ", "{:.4f}".format(q02 - q2))

        count = 0
        for i in range(corrections.size):
            if abs(corrections[i]) > conv_crit:
                cur = abs(corrections[i])
                count += 1
        convergence = count == 0
        itno += 1

    for i in range(busnum):
        if np.isnan(p_list[i]):
            p_list[i] = calcP(i, yBus, busnum, t_list, v_list)
        if np.isnan(q_list[i]):
            q_list[i] = calcQ(i, yBus, busnum, t_list, v_list)
    print('Final P and Q Values: ')
    for i in range(busnum):
        print("P", i + 1, ": ", "{:.4f}".format(p_list[i]), "\t\t\t", "Q", i + 1, ": ", "{:.4f}".format(q_list[i]))
    print('Final V and T Values: ')
    for i in range(busnum):
        print("V", i + 1, ": ", "{:.4f}".format(v_list[i]), "\t\t\t", "T", i + 1, ": ", "{:.4f}".format(t_list[i]))



#def distributionFactors(conv_crit, filename):
