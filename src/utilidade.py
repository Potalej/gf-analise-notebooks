from time import time
from src.arquivo import ler_arquivo_grande
import json
import matplotlib.pyplot as plt
from mecanica import mecanica
import numpy as np
import os
from scipy.stats import linregress

def infos_simulacao (pasta, chunksize=100)->dict:
    print(f"Lendo a pasta '{pasta}'... ", end='')

    # Lendo o arquivo "data.csv"
    t0 = time()
    arquivo_data = pasta + "/data.csv"
    qs, ps = ler_arquivo_grande(arquivo_data, chunksize)
    print(f'lido! ({round(time() - t0,2)}s)', end='')

    t0 = time()

    # Agora le o de valores iniciais
    with open(f"{pasta}/vi.json", 'r') as arq:
        json_obj = json.load(arq)

    # Algumas informacoes dos valores iniciais
    G = json_obj['G']
    N = json_obj['N']
    amortecedor = json_obj['integracao']['amortecedor']
    densidade = json_obj['colisoes']['densidade']
    massas = json_obj['valores_iniciais']['massas']
    M = sum(massas)
    checkpoints = json_obj['integracao']['checkpoints']
    q0 = np.array(json_obj['valores_iniciais']['posicoes'])
    p0 = np.array(json_obj['valores_iniciais']['momentos'])
    e0 = mecanica.energia_total(G, massas, q0, p0, amortecedor)

    # Agora vai calcular as informacoes para cada passo
    informacoes = {
        'intervalo': [json_obj['integracao']['t0'], json_obj['integracao']['tf'], json_obj['integracao']['checkpoints']],
        'i0': sum(massas[a]*(q0[a]@q0[a]) for a in range(N)),
        'e0': e0,
        'cinetica': [],
        'f_prod_q': [],
        'energia': [],
        'energia_erro': [],
        'virial': [],
        'virial_media': [],
        'vp': [],
        'vp_media': [],
        'autovalores': [[], [], []],
        'raio_meia_massa': [],
        'tempo_dinamico': [],
        'tempo_relaxacao_meia_massa': [],
        'centro_massas_minicentro': [],
        'raio_meia_massa_relativo': [],

        'escapes': [],
        'virial_central': [],
        'inercia_central': [],
        'complexidade': []
    }
    qs[0] = q0
    ps[0] = p0
    qs = np.array(qs)
    ps = np.array(ps)
    for i in range(checkpoints):
        # (Provaveis) escapes
        e_tots = mecanica.energia_total_separada(G, massas, qs[i], ps[i], amortecedor)
        positivos = 0
        indices = []
        for a, e_tot in enumerate(e_tots):
          if e_tot >= 0: positivos += 1
          else: indices.append(a)
        informacoes['escapes'].append(positivos/N)

        # Centro de massas
        # rcm = centro_massas(massas, np.array(qs[i]))
        # informacoes['centro_massas'].append(rcm)
        # p_total = np.zeros(3)
        # for p in ps[i]: p_total += p
        # informacoes['p_total'].append(np.linalg.norm(p_total))
        rcm_centro = mecanica.centro_massas(massas, np.array(qs[i]), indices)
        informacoes['centro_massas_minicentro'].append(np.linalg.norm(rcm_centro))

        # Com os que nao forem possiveis escapes, vamos considerar que sao o centro
        f_prod_qs = mecanica.virial_potencial_amortecido_separado(G,massas,qs[i],amortecedor)
        # momine_sep = mecanica.momento_inercia_separado(massas, qs[i])
        momine_central = mecanica.momento_inercia_indices(massas, qs[i], indices)
        virial_central = 0.0
        for a in indices:
          virial_central += (ps[i][a] @ ps[i][a])/massas[a] + f_prod_qs[a]
        informacoes['virial_central'].append(virial_central)
        informacoes['inercia_central'].append(momine_central)

        # Raio de meia massa relativo
        rmm_relativo = mecanica.raio_meia_massa(massas, np.array(qs[i]), rcm_centro)
        informacoes['raio_meia_massa_relativo'].append(rmm_relativo)

        # Raio de meia massa
        rmm = raio_meia_massa(massas, np.array(qs[i]))
        informacoes['raio_meia_massa'].append(rmm)

        # Tempo dinamico
        td = tempo_dinamico(rmm, G, M)
        informacoes['tempo_dinamico'].append(td)

        # Tempo de relaxacao da meia massa
        t_rh = tempo_relaxacao_meia_massa(len(massas), rmm, massas[0], G)
        informacoes['tempo_relaxacao_meia_massa'].append(t_rh)

        # Energia cinetica
        ec = mecanica.energia_cinetica(massas,ps[i])
        informacoes['cinetica'].append(ec)

        # Energia potencial
        f_prod_q, V = mecanica.termo_virial(G, massas, qs[i], amortecedor)
        informacoes['f_prod_q'].append(f_prod_q)

        # Energia total
        E = ec + V
        informacoes['energia'].append(E)
        informacoes['energia_erro'].append((E - e0))

        # Virial
        vir = ec + ec + f_prod_q
        informacoes['virial'].append(vir)

        # Virial (media)
        if i > 0: vir = (i*informacoes['virial_media'][-1] + vir)/(i+1)
        informacoes['virial_media'].append(vir)

        # Virial com potencial
        vir = ec + ec + V
        informacoes['vp'].append(vir)

        # Virial com potencial (media)
        if i > 0: vir = (i*informacoes['vp_media'][-1] + vir)/(i+1)
        informacoes['vp_media'].append(vir)

        # Tensor de inercia
        tensor = mecanica.tensor_inercia_geral(massas, qs[i])
        autovalores = np.linalg.eigvals(tensor)
        informacoes['autovalores'][0].append(autovalores[0])
        informacoes['autovalores'][1].append(autovalores[1])
        informacoes['autovalores'][2].append(autovalores[2])

        # Complexidade
        I = mecanica.momento_inercia(massas, qs[i])
        informacoes['complexidade'].append(-V*I)

    informacoes['cinetica'] = np.array(informacoes['cinetica'])
    informacoes['f_prod_q'] = np.array(informacoes['f_prod_q'])
    informacoes['energia'] = np.array(informacoes['energia'])
    informacoes['virial'] = np.array(informacoes['virial'])
    informacoes['virial_media'] = np.array(informacoes['virial_media'])
    informacoes['vp'] = np.array(informacoes['vp'])
    informacoes['vp_media'] = np.array(informacoes['vp_media'])
    informacoes['autovalores'] = np.array(informacoes['autovalores'])

    print(f' informacoes calculadas! ({round(time() - t0,2)}s)')

    return informacoes

def plot_info (titulo, x, y, dim=1, figsize=(8,2), log_scale_y=False):
    plt.figure(figsize=figsize)
    plt.title(titulo)

    if dim == 1: plt.plot(x,y)
    else:
        for yi in y: plt.plot(x,yi)

    if log_scale_y: plt.yscale('log')
    plt.show()

def eixo_x_intervalo (intervalo:list):
    return np.linspace(intervalo[0], intervalo[1], intervalo[2])

def plotar_informacoes (infos:dict, tf=0, figsize=(8,2))->None:
    intervalo = infos['intervalo']
    eixo_x = eixo_x_intervalo(intervalo)
    if tf > 0: condicao = (eixo_x <= tf)
    else: condicao = (eixo_x >= intervalo[0])

    # Plota conservacao da energia
    erro_energia = np.abs(infos['e0'] - infos['energia'])
    plot_info('Energia (erro)', eixo_x[condicao], erro_energia[condicao], figsize=figsize, log_scale_y=True)

    # Plota o virial (media)
    virial = infos['virial_media']
    plot_info('Virial (média)', eixo_x[condicao], virial[condicao], figsize=figsize)

    # Autovalores
    autovalores = [x[condicao] for x in infos['autovalores']]
    plot_info('Autovalores', eixo_x[condicao], autovalores, dim=3, figsize=figsize, log_scale_y=True)

def ler_json (arquivo:str):
    with open(arquivo, 'r') as arq:
        json_obj = json.load(arq)
    return json_obj

def ler_pastas_diretorio (diretorio:str)->dict:
    arquivos = dict()
    if diretorio[-1] != '/': diretorio = diretorio + '/'
    for pasta in os.listdir(diretorio):
        json_obj = ler_json(diretorio + pasta + '/vi.json')
        epsilon = json_obj['integracao']['amortecedor']
        if float(epsilon) == 0:
          densidade = json_obj['colisoes']['densidade']
          if json_obj['massas_iguais']:
            massa = json_obj['valores_iniciais']['massas'][0]
            densidade = round((3 * massa / (4 * np.pi * densidade))**(1/3), 4)

          contador = 0
          chave = str(densidade)
          while chave in list(arquivos.keys()):
            contador += 1
            chave = str(densidade) + "_" + str(contador)

          arquivos[chave] = {
              'pasta': diretorio + pasta,
              'json': json_obj
          }
        else:
          contador = 0
          chave = str(epsilon)
          while chave in list(arquivos.keys()):
            contador += 1
            chave = str(epsilon) + "_" + str(contador)
          arquivos[chave] = {
              'pasta': diretorio + pasta,
              'json': json_obj
          }
    return arquivos

def fazer_figura (arquivos:dict, chave:str, titulo:str, figsize:tuple=(8,3), lim_x=0, suavizar=False, y_referencia=[], plotar=True):
  plt.figure(figsize=figsize)
  plt.title(titulo)
  for epsilon in arquivos:
    infos = arquivos[epsilon]['infos']
    info = np.array(infos[chave])
    eixo_x = eixo_x_intervalo(infos['intervalo'])

    if lim_x == 0:
      eixo_x_plot = eixo_x
      eixo_y_plot = info
    else:
      eixo_x_plot = eixo_x[eixo_x <= lim_x]
      eixo_y_plot = info[eixo_x <= lim_x]

    if suavizar:
      eixo_y = np.zeros(len(eixo_y_plot))
      eixo_y[0] = eixo_y_plot[0]
      for i in range(1, len(eixo_y_plot)):
        eixo_y[i] = (i*eixo_y[i-1] + eixo_y_plot[i])/(i+1)
      eixo_y_plot = eixo_y

    plt.plot(eixo_x_plot, eixo_y_plot, label=epsilon)

  if len(y_referencia) > 0:
    for y_ref in y_referencia:
      plt.axhline(y_ref, c='black', linestyle='--')
  plt.legend()

  if plotar:
    plt.show()

def raio_meia_massa (massas, posicoes, denominador:float=2.0)->float:
  M = sum(massas)
  raios = [[i, q @ q] for i,q in enumerate(posicoes)]
  raios.sort()
  soma_massas = 0
  for i,r2 in raios:
    soma_massas += massas[i]
    if soma_massas >= M/denominador:
      return np.sqrt(r2)

def tempo_dinamico (rmm, G, M):
  rs = rmm/1.3
  return np.sqrt(rs / (G * M))

def tempo_relaxacao_meia_massa (N, rh, m, G=1, gamma=0.4):
  t_rh = 0.138 * np.sqrt((N*rh*rh*rh)/(G*m)) / np.log(gamma*N)
  return t_rh

def centro_massas (massas, posicoes):
  rcm=np.zeros(3)
  for m, q in zip(massas,posicoes): rcm += m*q
  return rcm/sum(massas)

def reta_mmq_local(x, y, xk, k=7):
    x = np.asarray(x)
    y = np.asarray(y)

    # 1. selecionar k vizinhos mais próximos
    idx = np.argsort(np.abs(x - xk))[:k]
    xi = x[idx]
    yi = y[idx]

    # 2. ajuste por MMQ (y = m x + b)
    slope, intercept, _, _, _ = linregress(xi, yi)
    return slope, intercept