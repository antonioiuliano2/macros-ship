#!/usr/bin/env python
# coding: utf-8

# # Test interfaccia per produrre tabella di vertici usando i file MC a disposizione
# Scrivo qui tutti i comandi in modo da non dimenticarmi quale script ho usato.
# La tabella di vertici è salvata in un file CSV (Comma Separated Values), con il primo rigo con gli indici
# Sunto degli script usati
# 1. Salvataggio indici di charm e figli di charm dalla simulazione MC
# 2. Identificazione vertici primari da MC
# 3. Identificazione vertici con figlie di charm
# 4. Studio extra-tracks
# 
# !**To do**: merge scripts 2 and 3 since they basically loop from the same file.

# In[1]:


#variables with file paths
macropath = "/home/antonio/Dottorato/Analisi/macros-ship"
workdir = "/home/antonio/Dottorato/Charmdata/CH1-R6/decay_search_MC"

#INPUT FILES
inputsimtree = workdir+"/ship.conical.Pythia8CharmOnly-TGeant4.root"
inputvertextree = workdir+"/vertextree_newformat.root"
inputtrackstree = workdir+"/linked_tracks.root"

#OUTPUT FILES
charmlistfile = workdir+"/charmlist.p"
outputcsv = workdir+"/MC_vertexlist_withallprimaries.csv"


# In[2]:


get_ipython().system('echo $charmlistfile')


# ## Salvataggio indici MC
# Per evitare di dover tenere aperto il file MC con la simulazione originaria, salvo delle liste con gli indici delle tracce di charm e figlie di charm in un file Python-Pickle. Salvo le seguenti liste:
# 
# 1. MCTrackID di adroni charmti
# 2. MCTrackID di figli di charm carichi
# 3. numero di figli di charm carichi
# 
# N.B. è fondamentale che gli altri script leggano queste liste **nello stesso ordine** in cui sono state scritte. Infatti Pickle non usa dei nomi per riconoscere gli oggetti salvati, ma li legge basandoli solo sull'ordine di scrittura e di lettura

# In[3]:


run $macropath"/analisi_charmsim/writecharmdaughters.py" -s $inputsimtree -co $charmlistfile


# ## Identificazione vertici primari
# Dato il formato della simulazione FairShip di ship-charm, le particelle prodotte al vertice primario hanno MotherID uguale a -1 (esclusi i due charm, che hanno MotherID 0). Pertanto posso riconoscere i primari nel MC richiedendo che la maggior parte delle tracce associate abbiamo MotherID -1. 
# 
# Se ci sono più vertici nello stesso evento con questa caratteristica (primari associati a vertici diversi), salvo solo il vertice con più tracce.

# In[4]:


run $macropath"/analisi_charmsim/search_primary_vertices.py" -f $inputvertextree -o $outputcsv


# ## Identificazione vertici secondari
# Ora utilizzo le liste salvate dallo script 1. Faccio un loop sulle tracce associate ai vertici: se MCTrackID corrisponde a una figlia di charm, salvo indici di traccia e vertice nel file csv in uscita.

# In[5]:


run $macropath"/analisi_charmsim/search_secondary_vertices.py" -f $inputvertextree -c $charmlistfile -o $outputcsv 


# ## Identificazione extra-tracks
# Non tutti i figli di charm sono associati a dei vertici. Alcuni sono collegati alla traccia madre (cioè l'adrone charmato), e quindi l'MCTrackID sarà riconosciuto come adrone charmato, essendo il primo segmento. Altri sono semplicemente non associati a nessun vertice. In entrambi i casi faccio un loop sul file di tracce:

# In[6]:


run $macropath"/analisi_charmsim/extratracks.py" -f $inputvertextree -t $inputtrackstree -c $charmlistfile -o $outputcsv 


# ## Controllo del file prodotto e ordinamento
# 
# Tutti gli script sono stati eseguiti. Ordino la lista per evento e topologia di vertice:
# 
# 1. Primari
# 2. Vertici con figlie di charm
# 3. Figli di charm collegati con il padre
# 4. Figli di charm non associati a nessun vertice
# 
# Posso controllare le prime e ultime righe per vedere se ci sono cose strane.

# In[7]:


import pandas as pd
df = pd.read_csv(outputcsv)


# In[8]:


df = df.sort_values(by=['MCEventID','topology'])


# In[9]:


df.to_csv(outputcsv,index=False)


# In[10]:


df


# In[ ]:


df[df['topology']==1] #primary


# In[ ]:


df[df['topology']==2] #secondary


# In[ ]:


df[df['topology']==3] #connected to charm parent


# In[ ]:


df[df['topology']==4] #not found in any vertex

