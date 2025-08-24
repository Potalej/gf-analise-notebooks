import requests
import zipfile

def baixar_zip (link:str, destino:str):
  response = requests.get(link)
  response.raise_for_status()
  with open(destino, 'wb') as f: f.write(response.content)
  print(f"Arquivo salvo com {destino}")

def descompactar (arquivo, destino:str):
  with zipfile.ZipFile(arquivo, "r") as zip_ref:
    zip_ref.extractall(destino)

  print("Extraído!")