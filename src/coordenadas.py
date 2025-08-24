import numpy as np

def parse_coordenadas(x, N):
    R = [
        [float(a) for a in x[0:N]],
        [float(a) for a in x[N:2*N]],
        [float(a) for a in x[2*N:3*N]]
    ]
    P = [
        [float(a) for a in x[3*N:4*N]],
        [float(a) for a in x[4*N:5*N]],
        [float(a) for a in x[5*N:6*N]]
    ]
    R = np.array(list(zip(*R)))
    P = np.array(list(zip(*P)))
    return R, P