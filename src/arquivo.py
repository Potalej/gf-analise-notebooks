import pandas as pd
from src.coordenadas import parse_coordenadas

def ler_arquivo_grande (arquivo, chunksize, saltos=1):
    valores_R, valores_P = [], []
    for chunk in pd.read_csv(arquivo, chunksize=chunksize, skiprows=3):
        # Aplica a funcao
        def func_res (x):
            [R, P] = [*parse_coordenadas(x, int(len(x)/6))]
            valores_R.append(R)
            valores_P.append(P)
        chunk.apply(func_res, axis=1)
    return valores_R[::saltos], valores_P[::saltos]