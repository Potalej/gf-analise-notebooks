from matplotlib.animation import FuncAnimation, FFMpegWriter
import matplotlib.pyplot as plt
from time import time
from src.arquivo import ler_arquivo_grande
import json
import numpy as np

def gerar_video_ffmpeg(eixo_t, pos, nome_arquivo="simulacao.mp4", fps=30, eixo_x=[], eixo_y=[], size=5):
    M, N, _ = pos.shape

    # Configuração inicial
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect("equal")

    # Determinar limites fixos
    if len(eixo_x) == 0:
      all_x = pos[:, :, 0].flatten()
      margem = 0.1 * max(all_x.max() - all_x.min())
      eixo_x = [all_x.min() - margem, all_x.max() + margem]

    if len(eixo_y) == 0:
      all_y = pos[:, :, 0].flatten()
      margem = 0.1 * max(all_y.max() - all_y.min())
      eixo_y = [all_y.min() - margem, all_y.max() + margem]

    ax.set_xlim(eixo_x)
    ax.set_ylim(eixo_y)

    scat = ax.scatter([], [], s=size, c="black")
    text = ax.text(eixo_x[0]+10, eixo_y[1]-10, '', ha='left', va='top', fontsize=12, color='black')

    # Função para atualizar cada quadro
    def update(frame):
        scat.set_offsets(pos[frame, :, :2])  # Apenas X e Y
        text.set_text(f"t = {eixo_t[frame]}")
        return scat, text

    ani = FuncAnimation(fig, update, frames=M, blit=True, interval=1000/fps)

    # Writer do ffmpeg (rápido e direto para MP4)
    writer = FFMpegWriter(fps=fps, codec="libx264", bitrate=-1)
    ani.save(nome_arquivo, writer=writer)

    plt.close(fig)

def gerar_video (pasta:str, nome_video:str, chunksize=100, fps=30, eixo_x=[], eixo_y=[], size=5)->dict:
    print(f"Lendo a pasta '{pasta}'... ", end='')

    # Lendo o arquivo "data.csv"
    t0 = time()
    arquivo_data = pasta + "/data.csv"
    qs, ps = ler_arquivo_grande(arquivo_data, chunksize)
    print(f'lido! ({round(time() - t0,2)}s)')


    # Lendo o json para gerar o eixo temporal
    with open(pasta + "/vi.json", 'r') as arq:
      jobj = json.load(arq)

    integracao = jobj['integracao']
    eixo_t = np.linspace(integracao['t0'], integracao['tf'], integracao['checkpoints'])

    t0 = time()
    gerar_video_ffmpeg(eixo_t, np.array(qs), f"{nome_video}.mp4", fps, eixo_x, eixo_y, size)
    print(f'vídeo {nome_video}.mp4 gerado! ({round(time() - t0,2)}s)')