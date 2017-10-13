import sys
import numpy as np
import math

def getX_K_Hat_Minus(x_k_hat_prev, u_k_prev):
    """
    Represents the time update equation for x_k_hat:
    x_k_hat_- = x_k-1_hat + u_k-1

    Parameters:
    x_k_hat_prev = x_k-1_hat
    u_k_prev = u_k-1

    Returns: x_k_hat as a 2-vector (2x1 numpy.matrix)
    """
    x_k_hat_minus = x_k_hat_prev + u_k_prev
    return x_k_hat_minus
    
def getP_K_Minus(P_k_prev):
    """
    Represents the time update equation for P_k_minus:
    P_k_minus = P_k-1 + Q

    Parameters:
    P_k 

    Returns: P_k_minus as a 2x2 matrix
    """
    Q11 = math.pow(10, -4) # For value at Q(1, 1)
    Q12 = 2 * math.pow(10, -5) # For value at Q(1, 2)
    Q = np.matrix("{0} {1}; {1} {0}".format(Q11, Q12))
    P_k_minus = P_k_prev + Q
    return P_k_minus
    
def getKalmanGain(P_k_minus):
    """
    Represents the equation for Kalman gain (K_k):
    K_k = (P_k_minus)/(P_k_minus + R)

    Parameters:
    P_k_minus = 2x2 matrix

    Returns: K_k as a 2x2 matrix
    """
    R11 = math.pow(10, -2) # For value at R(1, 1)
    R12 = 5 * math.pow(10, -3) # For value at R(1, 2)
    R = np.matrix("{0} {1}; {1} {0}".format(R11, R12))
    K_k = P_k_minus/(P_k_minus + R)
    return K_k

def getP_K(K_k, P_k_minus):
    """
    Represents the measurement update equation for P_k:
    P_k_minus = (I - K_k) * P_k_minus

    Parameters:
    K_k = 2x2 matrix
    P_k_minus = 2x2 matrix

    Returns: P_k as a 2x2 matrix
    """
    I = np.matrix("1 0; 0 1") # Identity matrix
    P_k = (I - K_k) * P_k_minus
    return P_k

def getX_K_Hat(x_k_hat_minus, K_k, z_k):
    """
    Represents the measurement update equation for x_k_hat:
    x_k_hat = x_k_hat_prev + K_k(z_k - x_k_hat_prev)

    Parameters:
    x_k_hat_minus = 2-vector
    K_k = 2x2 matrix
    z_k = 2-vector (given at beginning)

    Returns: x_k_hat, a 2-vector
    """
    x_k_hat = x_k_hat_minus + K_k * (z_k - x_k_hat_minus)
    return x_k_hat

if __name__ == "__main__":
    
    # Retrive file name for input data
    if(len(sys.argv) < 5):
        print "Four arguments required: python kalman2d.py [datafile] [x1] [x2] [lambda]"
        exit()
    
    filename = sys.argv[1]
    x10 = float(sys.argv[2])
    x20 = float(sys.argv[3])
    scaler = float(sys.argv[4])

    # Read data
    lines = [line.rstrip('\n') for line in open(filename)]
    data = []
    for line in range(0, len(lines)):
        data.append(map(float, lines[line].split(' ')))

    # Print out the data
    print "The input data points in the format of 'k [u1, u2, z1, z2]', are:"
    for it in range(0, len(data)):
        print str(it + 1) + " " + str(data[it])

    # Initial info
    I = np.matrix("1 0 ; 0 1")
    P_k_prev = scaler * I # Start with P_0
    x_k_hat_prev = np.matrix("{} ; {}".format(x10, x20)) # Start with x_0_hat

    averages_list = [] # List of x_k_hat 2-vectors found starting at k = 1
    variance_list = [] # List of P_k 2-vectors found starting at k = 1

    for k in range(len(data)):
        inputLst = data[k]
        u_1_k_prev = inputLst[0] # Same as u_1_k-1
        u_2_k_prev = inputLst[1] # Same as u_2_k-1
        u_k_prev = np.matrix("{} ; {}".format(u_1_k_prev, u_2_k_prev))
        z_1_k = inputLst[2]
        z_2_k = inputLst[3]
        z_k = np.matrix("{} ; {}".format(z_1_k, z_2_k))

        # Time update
        x_k_hat_minus = getX_K_Hat_Minus(x_k_hat_prev, u_k_prev)
        P_k_minus = getP_K_Minus(P_k_prev)

        # Kalman gain
        K_k = getKalmanGain(P_k_minus)

        # Measurement update
        P_k = getP_K(K_k, P_k_minus)
        x_k_hat = getX_K_Hat(x_k_hat_minus, K_k, z_k)

        # Store values found
        averages_list.append(x_k_hat)
        variance_list.append(P_k)

        # Setup variables for next iteration
        x_k_hat_prev = x_k_hat
        P_k_prev = P_k

    # Print the predicted x_k values and variances
    for k in range(len(averages_list)):
        print "x_{0}_hat = {1} \n P_{0} = {2}\n\n".format(k, averages_list[k], variance_list[k])
        
