import numpy as np
from scipy.integrate import tplquad

# 1. Definizione dell'integrando
# NOTA SciPy: tplquad richiede che l'ordine dei parametri sia f(z, y, x) -> f(phi, theta, p)
def integrando(phi, theta, p):
    # Evita problemi di divisione per zero se theta = pi/2
    sec_theta = 1.0 / np.cos(theta)
    
    # Flusso di muoni (o la tua funzione base)
    flusso = 18.0 / (p * np.cos(theta) + 145.0) * ((p + 2.7 * sec_theta)**(-2.7)) * ((p + 5.0) * np.sin(theta) / (p + 5.0 * sec_theta))
    
    # Fattore d'area / orientazione dipendente da phi e theta
    # Modifica questa riga se in futuro cambierà la dipendenza da phi!
    fattore_orientazione = np.abs(np.sin(phi) * np.sin(theta)) * np.sin(theta)
    
    return flusso * fattore_orientazione

# 2. Definizione dei limiti di integrazione
# p: da 1 a 10^5
p_min, p_max = 1.0, 1e5

# theta: da 0 a pi/2
# (usiamo un valore piccolissimo prima di pi/2 per evitare la singolarità di sec(theta))
theta_min = 0.0
theta_max = np.pi / 2.0 - 1e-9  

# phi: da 0 a 2*pi
phi_min = 0.0
phi_max = 2.0 * np.pi

# 3. Calcolo dell'integrale triplo
# Sintassi: tplquad(func, p_min, p_max, theta_min_func, theta_max_func, phi_min_func, phi_max_func)
risultato, errore_stimato = tplquad(
    integrando,
    p_min, p_max,                         # Limiti per p (variabile esterna x)
    lambda p: theta_min,                  # Limite inf per theta (variabile y)
    lambda p: theta_max,                  # Limite sup per theta (variabile y)
    lambda p, theta: phi_min,             # Limite inf per phi (variabile interna z)
    lambda p, theta: phi_max              # Limite sup per phi (variabile interna z)
)

print(f"Risultato Integrale Triplo: {risultato:.8f}")
print(f"Errore stimato da SciPy:   {errore_stimato:.2e}")