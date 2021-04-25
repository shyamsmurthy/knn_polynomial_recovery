# Creator Srinivas Vivek - January 2019
# Last Modified by Shyam S.M - April 2021
#
#
#@parallel(ncpus=2,timeout=150)
def divisorslist(p):
        return divisors(p)

# currently for 2^26 - takes about 1.5 mins or so
primelist=[]
def setup_prime_cache(x):
    p=2
    while p < 2^x:
        primelist.append(p)
        p=next_prime(p)


def divisorslist_bound(p,d,a,b):
    return reduce(list.__add__, ([i] for i in range(1, 2^b) if p % i == 0))

def old_divisorslist_bound_optimized(n,d,a,b): # NOTE to be used for less than 2^28 (b=28) only
    fl=[]
    bnd=28
    if b <= 28:
        bnd = b
    else:
        print "ERROR: use ecm factoring for beta > 28"
        return [1]

    if type(n) is not Integer:
        ny=Integer(int(n))
        yyf=ny.factor(limit=2^bnd,proof=False)
    else:
        yyf=n.factor(limit=2^bnd,proof=False)
    for i in range(len(yyf)):
        for j in range(yyf[i][1]):
            if yyf[i][0] <= 2^bnd:
                fl.append(yyf[i][0])

    flp = [[1]]
    row = 0
    ele = 0
    for i in range(0,len(fl)):
        if (fl[i] != ele):
            ele = fl[i]
            flp.append([ele])
            row += 1
            count = 1
        else :
            flp[row].append(ele ^ (count+1))
            count += 1

    div = [1]
    for i in range(1,len(flp)):
        tempdiv = []
        for j in range(0,len(flp[i])):
            for k in range(0,len(div)):
                if (div[k] * flp[i][j] < 2^b):
                    tempdiv.append(div[k] * flp[i][j])
        for j in range(0,len(tempdiv)):
            div.append(tempdiv[j])

    div = list(set(div))
    div.sort()
    return div


@fork(timeout=60,verbose=False)
def call_ecm_factor(n,b):
    EE=ECM()
    return EE.find_factor(n,B1=b)


import psutil
def kill_ecm_process():
    for process in psutil.process_iter():
        if process.pid == 0:
            continue
        if process.name() == 'ecm':
            try:
                process.kill()
            except psutil.NoSuchProcess :
                break
    return


def use_ecm_factoring(n,d,a,b):
    vbose = 0
    div = [1] # divisor list
    count = 0
    fl = [] # factor list
    copyn = n
    bnd=28

    if type(n) is not Integer:
#        print "(strange!!) type:", type(n), "for", n
        n=Integer(int(n))
        yyf=n.factor(limit=2^bnd,proof=False)
    else:
        yyf=n.factor(limit=2^bnd,proof=False)

    for i in range(len(yyf)):
        for j in range(yyf[i][1]):
            if yyf[i][0] <= 2^bnd:
                fl.append(yyf[i][0])
                n = n/yyf[i][0]

#    bnd_digits=Integer(2^bnd).ndigits()
    bnd_digits=100000
    while n > 2^bnd and n not in Primes():
        try:
            x=call_ecm_factor(n,bnd_digits)

        except OSError:
            print "OSError"
            continue # retry for same 'n'

        except AlarmInterrupt:
            print "AlarmInt"
            continue # retry for same 'n'

        if x[0] == "I":
            print "Factor x[0] : ", x[0]
            continue # retry for same 'n'

        if x[0] == "N": # we have timedout
            print "ECM Timeout, largest factor so far: ", fl[-1]
            kill_ecm_process()
            break

#        bnd_digits=Integer(x[0]).ndigits()
        if x[0] not in Primes():
            # Normally we should never get here !!
            print "Factor x[0] not prime!!"
            F = list(factor(x[0]))
            for i in range(0,len(F)):
                for j in range(0,F[i][1]):
                    if F[i][0] <= 2^b:
                        fl.append(F[i][0])
        else :
            if x[0] < 2^b:
                fl.append(x[0])

        n=n/x[0]
        if n != x[1]:
            # This should never happen!!
            print "Factors are different"

    # Finally, the residue if any...
    if n <= 2^b:
        print "We should have removed factors by now...", n, copyn
        if n in Primes():
            fl.append(n)

    fl.sort()
    if (vbose):
        print "Sorted fl:", fl

    flp = [[1]]
    row = 0
    ele = 0
    for i in range(0,len(fl)):
        if (fl[i] != ele):
            ele = fl[i]
            flp.append([ele])
            row += 1
            count = 1
        else :
            flp[row].append(ele ^ (count+1))
            count += 1

    for i in range(1,len(flp)):
        tempdiv = []
        for j in range(0,len(flp[i])):
            for k in range(0,len(div)):
                t = div[k] * flp[i][j]
                if t <= 2^b:
                    tempdiv.append(t)
        for j in range(0,len(tempdiv)):
            div.append(tempdiv[j])

    div = list(set(div))
    div.sort()
    return div



def get_consistency_array(finalstate,d,rn,stateF,skipF,vbose):

    consistent_array = [[0 for j in range(rn+1)] for i in range(2)]
    consistent_count = 0
    for i in range(0,rn+1):
        if finalstate[i] == 0:
            consistent_array[0][i] = 0
        else:
            consistent_array[0][i] = stateF[i][0][finalstate[i]]

    consistent_count = 1 # we will always send one set of divisors in [0]
    return consistent_array, consistent_count, 1



def get_valid_divisor_count(degree,rn,F,vbose):
    stateF = F
    numF = [0 for j in range(rn+1)] #num of non-zero factors
    skipF = [0 for j in range(rn+1)] #any row that needs to be skipped
    skipCount = 0
    for i in range(1,rn+1):
        for j in range(1,len(stateF[i][0])):
            if stateF[i][0][j] != 0:
                numF[i] += 1
            if numF[i] > 1:
                skipF[i] = 1
                skipCount += 1
    for i in range(1,rn+1):
        if numF[i] == 0:
            return -1,0,numF,skipF
    return 0,skipCount,numF,skipF


def compute_C_array_from_consistency_array(degree,rn,F,vbose,id):
    stateF = [[0 for j in range(rn+1)] for i in range(rn+1)]
    stateF = F
    C = [0]

    pivot_elem = stateF[1][0][id]
    C.append(pivot_elem)
    for i in range(2,degree+1):
        for j in range(1,len(stateF[i][0])):
            if stateF[i][0][j] == 0:
                continue
            tmp = stateF[i][0][j] - pivot_elem
            if tmp != 0 and tmp in stateF[i][1][0:len(stateF[i][1])]:
                C.append(stateF[i][0][j])
    return C


# Compute the F array factors that are consistent with one another
def compute_consistency_check_state(degree,rn,F,vbose):
    import time
    start_check_time = time.time()
    stateF = [[j for j in range(rn+1)] for i in range(rn+1)]
    non_unique = 0

    stateF = F
    finalstate = [0 for i in range(rn+1)]

    # Consistency of the set has already been verified - hence just clean up stateF and return consistent set

    flagc,skipCount,numF,skipF = get_valid_divisor_count(degree,rn,F,vbose)
    print "skipCount:", skipCount, "numF array:", numF, "skipF:", skipF
    for i in range(1,rn+1):
        if numF[i] > 1:
            non_unique = 1
            break
        for j in range(1,len(stateF[i][0])):
            if stateF[i][0][j] == 0:
                continue
            finalstate[i] = j

    if non_unique == 1:
        print "More than one set of divisors, will need deep consistency check..."
        return 0, -1, 0, 0, skipF

    if (flagc == -1 or (rn-skipCount) < degree+1):
        print "Unable to find d+1 divisors... rerun required.", degree, rn, skipCount
        return 0,-1, 0, 0, skipF

    consistent_array, consistent_count, flagc = get_consistency_array(finalstate,degree,rn,stateF,skipF,vbose)

    if flagc == -1:
        print "Error: No consistent set of divisors found..."
        return finalstate, -1, consistent_array, consistent_count,skipF

    return finalstate, 0, consistent_array, consistent_count,skipF



def internal_process_state_matrix(degree,col_id,rn,F,vbose,call_flag):
    stateF = [[j for j in range(rn+1)] for i in range(rn+1)]

    stateF = F
    vbose = 0
    start_row = col_id + 1 # use col id to get to right row
    next_row = start_row + 1

    if 1:
        rowFlag = 0
        for i in range(next_row,rn+1):  # index 0 is ignored
            nz_in = 0 # index 0 in every row has trivial divisor (=1) or 0

            while stateF[i-1][col_id][nz_in] == 0 or stateF[i-1][col_id][nz_in] == 1:
                nz_in += 1
                if nz_in == len(stateF[i-1][col_id]):
                    print "Encountered all 0 entry, problem row: ", i-1
                    print "Col id to examine", col_id
            least_nz_div = stateF[i-1][col_id][nz_in]

            # divisors in row i, that are less than the smallest non-zero div of previous row (i-1), are cleared
            if least_nz_div == 1:
                stateF[i][col_id][0] = 0
            cl_in = 0
            while stateF[i][col_id][cl_in] < least_nz_div:
                stateF[i][col_id][cl_in] = 0;
                cl_in += 1

            if cl_in != 0:
                 print "Row ", i, "no of elements zeroed : cl_in ", cl_in

            copy_nz  = least_nz_div

            for j in range(start_row,i):
                q = 1
                p = 1

                while (p < len(stateF[i][col_id])):
                    fFlag = 0
                    while (q < len(stateF[j][col_id])):
                        if ((stateF[i][col_id][p] != 0) and (stateF[j][col_id][q] != 0)):
                            rowFlag = 1
                            tmp = stateF[i][col_id][p] - stateF[j][col_id][q]
#                           if tmp != 0 and tmp in stateF[i][j][0:len(stateF[i][j])]:
                            if (tmp != 0) and (tmp in stateF[i][j]):
                                fFlag = 1

                        q += 1

                    if ((fFlag == 0) and (q == len(stateF[j][col_id]))): # we know stateF[i][col_id][p] is not compatible

                        stateF[i][col_id][p] = 0
                    q = 1
                    p += 1
                    while ((p < len(stateF[i][col_id])) and (stateF[i][col_id][p] == 0)) :
                        p += 1


                if (rowFlag == 0):
                    print call_flag,": Check fail row:", i, j,stateF[i][col_id],stateF[j][col_id]

            count_xyz = 0
            for xyz_abc in range(len(stateF[i][col_id])):
                if stateF[i][col_id][xyz_abc] != 0:
                    count_xyz += 1
            if count_xyz == 0:
                print "Problem state - ALL ZERO, copy_nz ", i, col_id, copy_nz
                print "p:", i, stateF[i][col_id]
                print "q:", j, stateF[j][col_id]
                print "in check:", stateF[i][j]

        for q in range(1,len(stateF[start_row][col_id])):
            if stateF[start_row][col_id][q] == 0:
                continue
            for i in range(next_row,rn+1):
                flag = 0
                rowNull = True
                for p in range(1,len(stateF[i][col_id])):
                    if (stateF[i][col_id][p] != 0):
                        rowNull = False # atleast one non-zero element in row
                        tmp = stateF[i][col_id][p] - stateF[start_row][col_id][q]
                        if tmp != 0 and tmp in stateF[i][col_id+1][0:len(stateF[i][col_id+1])] :
                            flag = 1

                if (rowNull == False) and (flag == 0): # this means qth element in first row is inconsistent
                    stateF[start_row][col_id][q] = 0

        rowFlag = 0
        for i in range(next_row,rn+1):  # index 0 is ignored
            for j in range(start_row,i):
                q = 1
                p = 1
                while (p < len(stateF[i][col_id])):
                    fFlag = 0
                    while (q < len(stateF[j][col_id])):
                        if ((stateF[i][col_id][p] != 0) and (stateF[j][col_id][q] != 0)):
                            rowFlag = 1
                            tmp = stateF[i][col_id][p] - stateF[j][col_id][q]
                            if tmp != 0 and tmp in stateF[i][j][0:len(stateF[i][j])]:
                                fFlag = 1

                        q += 1
                    if (fFlag == 0 and q == len(stateF[j][col_id])): # we know stateF[i][col_id][p] is not compatible
                        stateF[i][col_id][p] = 0
                    q = 1
                    p += 1
                    while ((p < len(stateF[i][col_id])) and (stateF[i][col_id][p] == 0)) :
                        p += 1
                if (rowFlag == 0):
                    print call_flag,": Check fail row:", i, j,stateF[i][col_id],stateF[j][col_id]

    return



# Check if we have the required number of consistent entries (degre+1)
def verify_consistency_check_count(degree,rn,F,vbose):
        stateF = [[j for j in range(rn+1)] for i in range(rn+1)]

        stateF = F
        vbose = 0
        if (vbose):
            print "Verify compute state Before: i:", 1, stateF[1][0]
            for i in range(2,rn+1):
                print "Verify compute state Before: i:", i, stateF[i][0], stateF[i][1]

        # clean up divisor list and try to find the unique and consistent divisor
        internal_process_state_matrix(degree,0,rn,F,vbose,22) # column 0

        numF = [0 for j in range(rn+1)] #num of non-zero factors
        skipF = [0 for j in range(rn+1)] #any row that needs to be skipped
        skipCount = 0
        for i in range(1,rn+1):
            for j in range(1,len(stateF[i][0])):
                if stateF[i][0][j] != 0:
                    numF[i] += 1
                if numF[i] > 3:
                    skipF[i] = 1
                    skipCount += 1
                    break
            if numF[i] == 0:
                skipF[i] = 1
                skipCount += 1

        if (vbose):
            print "skipCount:", skipCount, "numF array:", numF, "skipF:", skipF
            print "Verify compute state AFTER:i", 1, stateF[1][0]
            for i in range(2,rn+1):
                print "Verify compute state AFTER: i", i, stateF[i][0], stateF[i][1]

        if ((rn-skipCount) < degree+1):
            return -1, (rn-skipCount) # not found the needed number of entries
        else:
            return 0, (rn-skipCount)



# To check and make the matrix consistent for the given column and its adjacent (right) neighbor
def column_consistency_check(degree,col_id,rn,F,vbose):
        stateF = [[j for j in range(rn+1)] for i in range(rn+1)]

        stateF = F
        vbose = 0
        start_row = col_id + 1 # use col id to get to right row
        next_row = start_row + 1

        if (vbose):
            print "Col: Verify compute state Before: i:", start_row, stateF[start_row][col_id]
            for i in range(next_row,rn+1):
                print "Col: Verify compute state Before: i:", i, stateF[i][col_id], stateF[i][col_id+1]

        internal_process_state_matrix(degree,col_id,rn,F,vbose,55)

        if (vbose):
            print "Col: Verify compute state AFTER: i", start_row, stateF[start_row][col_id]
            for i in range(next_row,rn+1):
                print "Col: Verify compute state AFTER: i", i, stateF[i][col_id], stateF[i][col_id+1]

        return 0



def compute_deep_consistency_check(degree,a,b,rn,F,D,DCopy,vbose,index,column):
    print "Compute deep consistency check, column:", column
    stateF = F
    vbose = 0
    newF = [[] for i in range(rn+1)]

    for col_id in range(0,degree): # till degree
        row_id = col_id + 1    # first valid row in column 0 is 1 and so on : row_id = degree => col_id = degree - 1
        flist = [0 for i in range(row_id)] # all elems before first row are 0

        for i in range(row_id,degree+1): #degree + 1
            if vbose:
                print DCopy[i][col_id], "i, stateF", i, stateF[i]
            idx = 0
            for j in range(1,len(stateF[i][col_id])):
                if stateF[i][col_id][j] == 0:
                    continue
                if (idx == index):
                    flist.append(stateF[i][col_id][j])
                idx += 1
        for i in range (row_id,len(flist)):
            if flist[i] == 0:
                print "Divisor 0 in row", i, "in deep_consistency_check"
                return -1
        j = 1
        for i in range(row_id,degree+1):
            if vbose:
                print i, DCopy[i][col_id]/flist[i], flist[i], DCopy[i][col_id]
            DCopy[i][col_id] = DCopy[i][col_id]/flist[i]

        k = 1
        for i in range(row_id+1,degree+1):
            diff = DCopy[i][col_id] - DCopy[row_id][col_id]
            if (diff != 0):
                if diff <= 2^26:
                    tmp = divisorslist_bound_optimized(diff,degree,a,b)
                else:
                    tmp = use_ecm_factoring(diff,degree,a,b)
                SetTmp = Set(tmp)
                SetStateF = Set(stateF[i][col_id+1])
                SetInt = SetStateF.intersection(SetTmp)
                newList = SetInt.list()
                newList.append(1)
                newList.sort()
                if vbose:
                    print "",
                stateF[i][col_id+1] = newList
                DCopy[i][col_id+1] = diff

        for column_id in range(col_id+1,-1,-1):
            column_consistency_check(degree,column_id,rn,F,vbose)
            flagc,skipCount,numF, skipF = get_valid_divisor_count(degree,rn,F,vbose)
            if flagc == -1:
                print "RETURNING at get_valid_divisor_count...", flagc, skipCount, numF, skipF
                return -1

            if ((rn-skipCount) >= degree+1):
                print "",
                return 0

    return -1


# A subroutine to determine whether the given polynomial has a positive intege root

def check_posint_roots(tpoly):
        try:
                tl = tpoly.roots()
        except RuntimeError:
                return 0, -1
        for i in range(len(tl)):
                if ((tl[i][0].is_integer() == 1) and (tl[i][0] >= 0)):
                        return 1, tl[i][0]
        return 0, -1


def check_coefficient_bounds(tpoly):
## Now the only check is for integer coefficients
    abc = tpoly.coefficients(x0,sparse=False)
    for i in range(0,tpoly.degree(x0)+1):
        f = SR(abc[i])
        if (f.is_integer() != True):
            print "is_integer FALSE: ", f
            return -1, abc[i]
    return 0,1


# A subroutine to compute Lagrange polynomial based on possible divisor values
def compute_ai_coefficients_and_poly(d,x_diff,y_values,vbose):
    polyX = 0
    vbose = 0
    if (vbose):
        print len(x_diff), x_diff
        print "d:", d, "y_val len:", len(y_values), "y val:", y_values

    numerator = 1
    for i in range(0,d+1) :
        numerator *= (x - x0 - x_diff[i])

    print "Coefficients from constant to highest degree"
    L = 0
    for i in range(0,d+1):
        denom = 1
        for j in range(0,d+1) :
            if (i == j): continue
            denom = denom * (x_diff[i] - x_diff[j])
        L += (y_values[i]/denom)*(numerator/(x - x0 - x_diff[i]))

    tpoly = L.coefficients(x,sparse=False)
    print d
    for xyz in range(0,d+1):
        print tpoly[xyz].coefficients(x0,sparse=False), ","
        print " "

    for i in range(0,d+1):
        if tpoly[i].degree(x0) == 1:
            degree1poly = i
        flag,coeff = check_coefficient_bounds(tpoly[i])

        if flag == -1:
            if (vbose): print "Discarding:", x_diff
            return -1, polyX

    x0_val = tpoly[degree1poly].roots(x=x0) # compute the root of linear polynomial
    intx0_val = floor(x0_val[0][0])

    # create list of coefficients evaluated at x0
    polyX = 0
    a = []
    for i in range(0,d):
        a.append(tpoly[i](x0=intx0_val))
    a.append(tpoly[d])

    for i in range(0,d+1):
        polyX += a[i]*x^i
    print "Recovered: x0", intx0_val, "and Polynomial:", polyX

    return 0, polyX


var('x,x0')
def enum_over_all_consistency_elements(d,rn,b,state,F,P,X,vbose,consistency_array,consistency_count,skipF) :
    import time
    discard_count = 0
    start_enum_time = time.time()
    R = PolynomialRing(QQ, 'x')

    print "consistency_array:", consistency_array, " count: ", consistency_count, "skipF:", skipF

    y_output = []
    y_index = 0
    y_count = 0
    while y_count < d+1:
        if skipF[y_index] != 1:
            y_output.append(P[y_index])
            y_index += 1
            y_count += 1
        else:
            y_index += 1

    for i in range(0,d+1):
        if skipF[i] != 1:
            y_output.append(P[i])

    failflag = 0
    for i in range(0,consistency_count):
        fc, polyX = compute_ai_coefficients_and_poly(d,consistency_array,y_output,vbose)
        outP = []
        if (fc == -1):
            discard_count += 1
        if (polyX != 0):
            for j in range(0,len(X)):
                outP.append(polyX(x=X[j]))
            for j in range(0,len(P)):
                if P[j] != outP[j]:
                    failflag = 1
                    print "Outputs differ at:", j, "for",polyX
                    break # if one fails then stop

    if (failflag == 0) and (discard_count == 0):
        print "Successfully verified for all", len(X)+1, "inputs."
        print "Consistent divisor set got:", consistency_count, "discarded:", discard_count
    else:
        print "Unsuccessful: discard_count:", discard_count


def gen_matrix_polyeval(d,a,b,n,rn,bnd,vbose,given_poly,polyOpArray,ipArray):
        import time

        if given_poly != 0 : # This is the case where the polynomial outputs are given
                P = [(long)(polyOpArray[i]) for i in range (n)]
                A = [0 for i in range(d+1)]
                X = [ipArray[i] for i in range(n)]
        else : # otherwise we pick uniform random
                A = [randint(1,(2^a)-1) for i in range(d+1)]

                X = [randint(0,(2^b)-1) for j in range(n)]
                X = list(set(X)) # make the list unique
                X.sort()        # This is not a problem as the polynomial is monotonic
                n = len(X)
                P = [ sum(A[i]* X[j]^i for i in range(d+1)) for j in range(n) ] # vector of outputs

        D = [[0 for j in range(n)] for i in range(n)]
        DCopy = [[0 for j in range(n)] for i in range(n)]
        F = [[] for i in range(n)]
        nF = [[] for i in range(n)]

        print "X[0:rn]:", X[0:rn]
        gen_mx_start = time.time()

        dbg_val=0
        for i in range(1,rn+1):
                for j in range(0,i):
                        exception_count = 0
                        D[i][j] = P[i] - P[j]
                        DCopy[i][j] = D[i][j]
                        if(D[i][j]!=0):
                            print i,j,D[i][j]
                            vtime=time.time()
                            if (b <= 28): # for small bounds use vanilla factorisation
                                tmp = divisorslist_bound_optimized(D[i][j],d,a,b)
                            else:
                                tmp = use_ecm_factoring(D[i][j],d,a,b)
                            if (X[i] - X[j]) not in tmp:
                                print "FACTOR NOT PRESENT:", X[i] - X[j]

                            if dbg_val < 5:
                                dbg_val += 1
                                print time.time()-vtime
                            F[i].append(tmp)
                            nF[i].append(len(tmp))
                            if (j == 0) and (vbose == 1):
                                print "Div List[0] for i",i, len(F[i][0]), F[i][0]
                        else:
                            F[i].append([0])
                            nF[i].append(len([0]))

                if i == 1:
                    continue

                if (i >= d+2):
                    vbose = 0
                    column_consistency_check(d+1,1,i,F,vbose)
                    flagV,count = verify_consistency_check_count(d+1,i,F,vbose)
                    if flagV == -99:
                        exit()
                    if (count >= d+1):
                        print "Obtained needed x-differences:", count
                        break
                    else:
                        print "Iteration i:", i, "count:", count

        if (vbose):
            print "Exception count", exception_count, "Step 0 done. time: ", time.time() - gen_mx_start
            print "A:", A[::-1], "\n"
            print "X:", X, "\n"
            print "P:", P, "\n"

        return A,X,P,F,D,DCopy # return params to prevent error message from caller



def solve_poly_eval(d,a,b,n,rn,bnd,vbose,given_poly,polyOpArray,ipArray):

    import time
    import sys

    start_time = time.time()

# We are genrating A,X only for testing purposes.
    R = PolynomialRing(QQ, 'x')

    A,X,P,F,Diff,DiffCopy = gen_matrix_polyeval(d,a,b,n,rn,bnd,vbose,given_poly,polyOpArray,ipArray)
    gen_matrix_time = time.time()

    div_index = 0
    fail_count = 0
    state, flagc,consistency_array, consistency_count, skipF = compute_consistency_check_state(d+1,d+2,F,vbose)

    c_count = 0
    CC = []
    for i in range (1,len(F[1][0])):
        if F[1][0][i] != 0:
            c_count += 1
            C = compute_C_array_from_consistency_array(d,rn,F,vbose,i)
            CC.append(C)

    consistency_check_time = time.time()
    print "x0 chosen:", X[0], "Random Polynomial:", A[d:0:-1], A[0] # print in reverse
    print "A:", A
    print "Poly Gen time:", gen_matrix_time - start_time

    for ii in range(d+2):
        skipF[ii] = 0
    for kk in range(c_count):
        enum_time = time.time()
        print "For CC[]: ", kk
        enum_over_all_consistency_elements(d,rn,b,state,F,P,X,vbose,CC[kk],1,skipF);
        print "Time (poly gen + this iter):", gen_matrix_time - start_time + time.time() - enum_time
        print "######"

    return 1



def main(a,b,d):
    PGiven = PolynomialRing(QQ, 'x')
    PGiven = 0*x+0

    n = 20 # n gets changed below to be 3 * degree
    polyOpArray = [ 0 for i in range(n) ] # vector of outputs

    polyOpArraySort = []
    ipArraySort = []
    min_degree = 10
    max_degree = 12
    degree = max_degree
    iter = 1
    alpha = 128
    beta  = 64
    print "Start Execution"

    alpha=a
    beta=b
    degree=d
    n = 2 * degree
    rn = degree + 5 # factors will be computed for rn elements
    bnd = rn  # bound on number of y values to pick - set to rn for now
    vbose = 0

    print "alpha", a, "beta", b, "degree", d
    if 1:
        for i in range(0,iter):
            print "Starting Execution for degree", degree, "polynomial for (alpha, beta):", alpha, beta
            for j in range(1,iter+1):
                print "Run id :", j, "(",degree,alpha,beta,")"
                solve_poly_eval(degree,alpha,beta,n,rn,bnd,vbose,PGiven,polyOpArraySort,ipArraySort)
                print "\n"
            degree += 0 # was 50
            n = 2*degree
            rn = degree + 5
            bnd = rn
            beta += 0
            print "###########"

#
#
# To run give main(alpha, beta, degree)
#
# Example: main(128,32,32)
#
#
