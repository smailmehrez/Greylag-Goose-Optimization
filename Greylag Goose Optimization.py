from Fonction import*
import numpy as np
import random
import matplotlib.pyplot as plt
import math
import warnings
warnings.filterwarnings('ignore')

def GGO(func,Solution,nb_solu,taille_solu,itération):
    T = itération
    taille_solution = taille_solu
    nb_solution = nb_solu
    Matrice_Solution=np.copy(Solution)

    # ******** Trouver the best Solution and the worst *************
    idxMin = np.argmin(Matrice_Solution[:,taille_solution])
    solution_best_in_itra=np.copy(Matrice_Solution[idxMin,:])
    vecteur_Fonc_objectiv_meillure=np.full(shape=(T+1), fill_value=Matrice_Solution[idxMin,taille_solution],dtype=float)
    _,UB,_,LB=func(([1,2]))
    #******Inialiser les parametre utilise dans algorithme GGO HHO******
    #******************GGO******************************************
    t = 1
    n1 = int( nb_solution/2)
    n2 = int(nb_solution/2)
    a = np.linspace(0, 2, num=T)
    while t <= T:
        z = 1-(t/T)**2 
        r1 = random.uniform(0,1)
        r2 = random.uniform(0,1)
        r3 = random.uniform(0,1)
        r4 = random.uniform(0,1)
        r5 = random.uniform(0,1)
        A = 2 * a[t-1] * r1 - a[t-1]
        C = 2 * r2
        w1 = random.uniform(0,2)
        w2 = random.uniform(0,2)
        w3 = random.uniform(0,2)
        w4 = random.uniform(0,2)
        w = random.uniform(0,2)
        b = 10
        l = random.uniform(-1,1)
        for i in range(n1):
            if (t % 2) == 0 :
                if r3 < 0.5 :
                    if abs(A) < 1:
                        Matrice_Solution[i,:-1] = solution_best_in_itra[:-1] -A * abs(C*solution_best_in_itra[:-1]-Matrice_Solution[i,:-1])
                    else :
                        Sa = random.randint(0,n1-1)
                        Sb = random.randint(0,n1-1)
                        Sc = random.randint(0,n1-1)
                        while (Sa == idxMin) or (Sb == idxMin) or (Sc == idxMin) :
                            Sa = random.randint(0,n1-1)
                            Sb = random.randint(0,n1-1)
                            Sc = random.randint(0,n1-1)
                            while Sb == Sc :
                                Sb = random.randint(0,n1-1)
                        Matrice_Solution[i,:-1] = w1 * Matrice_Solution[Sa,:-1] + z * w2 * (Matrice_Solution[Sb,:-1] - Matrice_Solution[Sc,:-1]) + (1 - z) * w3 * (Matrice_Solution[i,:-1] - Matrice_Solution[Sa,:-1])
                else :
                    Matrice_Solution[i,:-1] = w4 * abs(solution_best_in_itra[:-1] - Matrice_Solution[i,:-1]) * math.exp(b*l) * math.cos(2*math.pi*l) + (2*w1*(r4+r5)) * solution_best_in_itra[:-1]
            else :
                Matrice_Solution[i,:-1] = Matrice_Solution[i,:-1] + (1+z) * w *(Matrice_Solution[i,:-1] - solution_best_in_itra[:-1])
            for j in range(taille_solution):
                if Matrice_Solution[i,j] < LB:
                    Matrice_Solution[i,j] = "{0:.3f}".format(random.uniform(LB,0))
                if Matrice_Solution[i,j] > UB:
                    Matrice_Solution[i,j] = "{0:.3f}".format(random.uniform(0,UB))
            Matrice_Solution[i,taille_solution],_,_,_=func(Matrice_Solution[i,:-1])
        for i in range(n2):
            if (t % 2) == 0 :
                x1 = random.randint(0,n2-1)
                x2 = random.randint(0,n2-1)
                x3 = random.randint(0,n2-1)
                S1 = Matrice_Solution[n1+x1,:-1] - A * abs(C * Matrice_Solution[n1+x1,:-1] - Matrice_Solution[n1+i,:-1])
                S2 = Matrice_Solution[n1+x2,:-1] - A * abs(C * Matrice_Solution[n1+x2,:-1] - Matrice_Solution[n1+i,:-1])
                S3 = Matrice_Solution[n1+x3,:-1] - A * abs(C * Matrice_Solution[n1+x3,:-1] - Matrice_Solution[n1+i,:-1])
                Matrice_Solution[n1+i,:-1] = (S1+S2+S3)/3
            else:
                Matrice_Solution[n1+i,:-1] = Matrice_Solution[n1+i,:-1] + (1+z) * w *(Matrice_Solution[n1+i,:-1] - solution_best_in_itra[:-1])
            for j in range(taille_solution):
                if Matrice_Solution[n1+i,j] < LB:
                    Matrice_Solution[n1+i,j] = "{0:.3f}".format(random.uniform(LB,0))
                if Matrice_Solution[n1+i,j] > UB:
                    Matrice_Solution[n1+i,j] = "{0:.3f}".format(random.uniform(0,UB))
            Matrice_Solution[n1+i,taille_solution],_,_,_=func(Matrice_Solution[n1+i,:-1])

        idxMin = np.argmin(Matrice_Solution[:,taille_solution])
        if (solution_best_in_itra[taille_solution]> Matrice_Solution[idxMin,taille_solution]):
            solution_best_in_itra=np.copy(Matrice_Solution[idxMin,:])
        vecteur_Fonc_objectiv_meillure[t] = solution_best_in_itra[taille_solution]

        if (vecteur_Fonc_objectiv_meillure[t] >= vecteur_Fonc_objectiv_meillure[t-1] or vecteur_Fonc_objectiv_meillure[t] >= vecteur_Fonc_objectiv_meillure[t-2]) and n1 < nb_solution-1:
            n1=n1+1
            n2=n2-1
        else:
            n2=n2+1
            n1=n1-1 
        t = t + 1
    return(vecteur_Fonc_objectiv_meillure)
                


T=100
n=10
fonction=Quartic
_,UB,taille_solution,LB=fonction(([1,2]))
nb_solution=n
# crer matrice du Solution
M1 = np.full(shape=(nb_solution, taille_solution+1), fill_value=0,dtype=np.float64)
# inialiser le Matrice du solution random
for j in range(nb_solution):
    for i in range(taille_solution):
        M1[j,i] = "{0:.3f}".format(random.uniform(LB, UB)) 
# calcule la fonction objective 
for j in range(nb_solution):
    M1[j,taille_solution],_,_,_=fonction(M1[j,:-1])
vecteur_best=GGO(fonction,M1,nb_solution,taille_solution,T)
ax1 =plt.axes()
ax1.plot(vecteur_best,label="GGO",c="#4B2991")
ax1.set_ylabel('Coût F')
ax1.set_xlabel("Nombre d'itération")
if fonction== Ackley or fonction == Michalewicz or fonction== Schwefel1 or fonction== Schwefel2_26 or fonction== Himmelblau:
    ax1.legend()
    ax1.set_title("Le Coût f en fonction de Nombre d'itération")
    plt.show()
else:
    ax1.set_yscale('log')
    ax1.legend()
    ax1.set_title("Le Coût f en fonction de Nombre d'itération")
    plt.show()