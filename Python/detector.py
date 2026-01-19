import numpy as np
import itertools
from mapping import qam_symbols_detection

def zf_detector(H, y, constellation):
    '''
    Input
        - H: channel matrix with shape (Nr, Nt)
        - y: received symbol vector with shape (Nr, 1)
        - constellation: constellation used for transmitted symbols
    Output
        Zero Forcing (ZF) detected symbols
    '''
    W = np.linalg.inv(H.conj().T @ H) @ H.conj().T
    x_hat = W @ y 
    detected_symbols = qam_symbols_detection(symbols=x_hat.flatten(), constellation=constellation)

    return detected_symbols

def mmse_detector(H, y, SNR, constellation):
    '''
    Input
        - H: channel matrix with shape (Nr, Nt)
        - y: received symbol vector with shape (Nr, 1)
        - SNR: linear SNR
        - constellation: constellation used for transmitted symbols
    Output
        Minimum Mean Square Error (MMSE) detected symbols
    '''

    Nt = H.shape[1]
    W = np.linalg.inv(H.conj().T @ H + np.identity(Nt) / SNR) @ H.conj().T
    x_hat = W @ y
    detected_symbols = qam_symbols_detection(symbols=x_hat.flatten(), constellation=constellation)

    return detected_symbols

def maximum_likelihood_detector(H, y, constellation):
    '''
    Input
        - H: channel matrix with shape (Nr, Nt)
        - y: received symbol vector with shape (Nr, 1)
        - constellation: constellation used for transmitted symbols
    Output
        - x_hat: equalized transmitted symbol vector with shape (Nt, 1)
    '''
    Nt = H.shape[1]

    # Generate M ** Nt possible transmitted vector (M: modulation order of transmitted symbols)
    x = np.array(list(itertools.product(constellation, repeat=Nt))).T  # shape: (Nt, M ** Nt)
    y_ref = H @ x   # shape: (Nr, M ** Nt)

    # Compare the distance between y and each possible y_ref = H @ x, and find y_ref[:,idx] that is nearest to y
    dist = np.sum(np.abs(y - y_ref) ** 2, axis=0)  # shape: (M ** Nt,)
    idx = np.argmin(dist)

    x_hat = x[:, idx].reshape(Nt, 1)

    return x_hat

def simple_sic_detector(H, y, constellation):
    '''
    Perform simple Successive Interference Cancellation (SIC) equalizer with random decoding order.
    'sic_equalizer' will detect symbol in strong-to-weak channel order.

    Input
        - H: channel matrix with shape (Nr, Nt)
        - y: received symbol vector with shape (Nr, 1)
        - constellation: constellation used for transmitted symbols
    Output
        - x_hat: detected transmitted symbol vector with shape (Nt, 1)
    '''
    Nt = H.shape[1]

    # Q: (Nr, Nr), R: (Nr, Nt), H: (Nr, Nt)
    Q, R = np.linalg.qr(H, mode='complete')

    # Q^H @ y = R @ x + Q^H @ n
    Qy = Q.conj().T @ y

    x_hat = np.zeros(shape=(Nt, 1), dtype=np.complex128)

    for sym_i in range(Nt - 1, -1, -1):
        # Q^H @ y (LHS) = R @ x(RHS)
        rhs_sum = 0
        for i in range(sym_i, Nt - 1):
            rhs_sum += R[sym_i, i + 1] * x_hat[i + 1, 0]
        x_hat[sym_i, 0] = (Qy[sym_i, 0] - rhs_sum) / R[sym_i, sym_i]

        x_hat[sym_i, 0] = constellation[np.argmin(np.abs(constellation-x_hat[sym_i, 0]))]

    return x_hat

def sorted_QR(H):
    '''
    Given channel matrix H and provided sorted Q, R and permutation matrices

    Input
        - H: channel matrix with shape (Nr, Nt)
    Output
        - Q: unitary matrix after sorted QR decomposition (Nr, Nt)
        - R: upper triangular matrix after sorted QR decomposition (Nt, Nt)
        - P: permutation matrix that sorts H in weak-to-strong (column) order (Nt, Nt)
        - order: weak-to-strong order of the columns in matrix H (Nt,)
    '''
    Nr = H.shape[0]
    Nt = H.shape[1]

    # H @ P = Q @ R
    Q = H.copy().astype(np.complex128)
    R = np.zeros((Nt, Nt), dtype=np.complex128)
    P = np.eye(Nt, dtype=int)

    for i in range(Nt):
        s = np.real(np.diag(Q[:,i:].conj().T @ Q[:,i:]))   # s[i] is the channel sum of transmitter i (sum_j h_ji^2 )
        col_idx = np.argmin(s)
        col_idx = col_idx + i

        # swap column 'col_idx' to column i
        temp = Q[:,i].copy()
        Q[:,i] = Q[:, col_idx]
        Q[:, col_idx] = temp

        temp = P[:,i].copy()
        P[:,i] = P[:, col_idx]
        P[:, col_idx] = temp

        temp = R[:,i].copy()
        R[:,i] = R[:, col_idx]
        R[:, col_idx] = temp    

        R[i, i] = np.sqrt(s[col_idx-i])
        Q[:, i] = Q[:, i] / R[i, i]

        # interference cancellation: remove ith column's component in the remaining columns
        for k in range(i+1, Nt):
            R[i, k] = Q[:, i].conj().T @ Q[:, k]
            Q[:, k] = Q[:, k] - R[i, k] * Q[:, i]
    
    order = P.argmax(axis=0)

    return Q, R, P, order

def sic_detector(H, y, constellation):
    '''
    Perform simple Successive Interference Cancellation (SIC) equalizer which detect symbol in strong-to-weak channel order.

    Input
        - H: channel matrix with shape (Nr, Nt)
        - y: received symbol vector with shape (Nr, 1)
        - constellation: constellation used for transmitted symbols
    Output
        - x_hat: equalized transmitted symbol vector with shape (Nt, 1)
    '''

    Nt = H.shape[1]

    Q, R, P, order = sorted_QR(H=H)

    # y = H @ x + noise
    # y = QR @ P^H @ x + noise
    # Q^H @ y = R @ (P^H @ x) + Q^H @ n
    Qy = Q.conj().T @ y

    perm_x_hat = np.zeros(shape=(Nt, 1), dtype=np.complex128)  # P^H @ x

    for sym_i in range(Nt - 1, -1, -1):
        # Q^H @ y (LHS) = R @ x(RHS)
        rhs_sum = 0
        for i in range(sym_i, Nt - 1):
            rhs_sum += R[sym_i, i + 1] * perm_x_hat[i + 1, 0]
        perm_x_hat[sym_i, 0] = (Qy[sym_i, 0] - rhs_sum) / R[sym_i, sym_i]
        perm_x_hat[sym_i, 0] = constellation[np.argmin(np.abs(constellation-perm_x_hat[sym_i, 0]))]

    x_hat = P @ perm_x_hat

    return x_hat

def mmse_sic_detector(H, y, SNR, constellation):
    '''
    Perform simple Successive Interference Cancellation (SIC) equalizer which detect symbol in strong-to-weak channel order.

    Input
        - H: channel matrix with shape (Nr, Nt)
        - y: received symbol vector with shape (Nr, 1)
        - SNR: linear SNR
        - constellation: constellation used for transmitted symbols
    Output
        - x_hat: equalized transmitted symbol vector with shape (Nt, 1)
    '''

    Nr = H.shape[0]
    Nt = H.shape[1]

    alpha = np.sqrt(Nt / SNR)
    H_aug = np.concatenate((H, alpha*np.eye(Nt)), axis=0)
    Q_aug, R_aug, P, order = sorted_QR(H=H_aug)
    Q = Q_aug[:Nr,:Nt]
    R = R_aug[:Nt,:]

    # y = H @ x + noise
    # y = QR @ P^H @ x + noise
    # Q^H @ y = R @ (P^H @ x) + Q^H @ n
    Qy = Q.conj().T @ y

    perm_x_hat = np.zeros(shape=(Nt, 1), dtype=np.complex128)  # P^H @ x

    for sym_i in range(Nt - 1, -1, -1):
        # Q^H @ y (LHS) = R @ x(RHS)
        rhs_sum = 0
        for i in range(sym_i, Nt - 1):
            rhs_sum += R[sym_i, i + 1] * perm_x_hat[i + 1, 0]
        perm_x_hat[sym_i, 0] = (Qy[sym_i, 0] - rhs_sum) / R[sym_i, sym_i]
        perm_x_hat[sym_i, 0] = constellation[np.argmin(np.abs(constellation-perm_x_hat[sym_i, 0]))]
    x_hat = P @ perm_x_hat

    return x_hat