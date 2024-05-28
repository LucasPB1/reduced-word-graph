#!/usr/bin/env sage

# author : Lucas POUILLART

import time

# On va commencer par redéfinir la fonction de réduction pour les groupes de Weyl, afin de l'appliquer
# aux groupes de Coxeter en général :

# On utilise le fait que dans sage, la fonction roots() renvoie une liste contenant d'abord les racines
# positives puis les négatives
def longueur(sigma, W) :
    """
    Renvoie la longueur de l'élément sigma dans le système de Coxeter W

    La longueur de l'élément se calcule simplement en voyant combien de
    racines positives sont envoyées sur des racines négatives par l'action de sigma.
    (Voir le lemme 2.8)

    EXEMPLE ::
    
        sage: W = CoxeterGroup(['A',3])
        sage: S = W.gens()
        sage: w = S[2]*S[1]*S[2]*S[2]*S[0]*S[2]*S[1]*S[2]
        sage: longueur(w,W)
        6
    """
    pi = W.positive_roots()
    minus_pi = W.roots()[len(pi):]
    res = [] #liste qui contiendra l'intersection de pi et de sigma ^-1 (-pi)
    #boucle qui teste pour tout alpha dans pi si sigma(alpha) dans minus_pi, et si oui l'ajoute à res
    for alpha in range(len(pi)) :
        beta = W.roots()[sigma.action_on_root_indices(alpha)]
        if beta in minus_pi :
            res += [alpha]
    return len(res)

def associate_root(sigma,W):
    """
    Pour une réflexion sigma appartenant à W, renvoie l'indice de la racine de W
    associée à la réflexion sigma

    Fonctionne pour des groupes de réflexions ou des groupes de Coxeter
    par leur représentation géométrique

    EXEMPLE ::
    
        sage: W = CoxeterGroup(['A',3])
        sage: S = W.gens()
        sage: w = S[2]
        sage: associate_root(w,W)
        5
        sage: W.roots()[5]
        (0, 0, 1)
    """
    n = len(W.positive_roots())
    for alpha in range(n):
        k = sigma.action_on_root_indices(alpha)
        if W.roots()[k] == - W.positive_roots()[alpha] :
            return k - n
    return -1

# On écrit une fonction permettant la construction des
# applications partielles sigma[i+1] ... sigma[j-1]
def constructPartialSigma(sigma, i, j, W):
    w = W.one()
    for s in sigma[i+1:j]:
        w = w * s
    return w

# On adapte la fonction de réduction du mot aux groupes de Coxeter en général
# ce qu'on justifie par le fait que tout groupe de Coxeter admet une représentation
# fidèle comme groupe de réflexions réel

def reduction(sigma,W) :
    """
    Pour une liste de réflexions simples sigma représentant un mot sur S, retourne
    un mot réduit équivalent par l'application successive du théorème de condition
    de délétion (voir théorème 2.7)

    Fonctionne pour les groupes de Coxeter finis (le critère de longueur théoriquement 
    vrai dans tous les groupes de Coxeter ne peut pas s'appliquer informatiquement 
    dans les groupes infinis)

    EXEMPLE ::
    
        sage: W = CoxeterGroup(['A',3])
        sage: S = W.gens()
        sage: w = [S[2],S[0],S[0],S[1],S[0],S[1],S[1],S[2],S[1],S[2]]
        sage: reduction(w,W)
        [
        [ 1  0  0]  [ 1  0  0]  [-1  1  0]  [ 1  0  0]  [ 1  0  0]  [ 1  0  0]
        [ 0  1  0]  [ 1 -1  1]  [ 0  1  0]  [ 0  1  0]  [ 1 -1  1]  [ 0  1  0]
        [ 0  1 -1], [ 0  0  1], [ 0  0  1], [ 0  1 -1], [ 0  0  1], [ 0  1 -1]
        ]
    """
    w = W.one()
    b = False
    sigma2 = sigma.copy()
    for s in sigma :
        w = w * s
    l = longueur(w,W)
    while len(sigma2) > l :
        j = 1
        while j < len(sigma2) and not(b):
            alpha = associate_root(sigma2[j],W)
            i = 0
            while i < j and not(b):
                w = constructPartialSigma(sigma2, i, j, W)
                alphaI = w.action_on_root_indices(alpha)
                if associate_root(sigma2[i],W) == alphaI:
                    del(sigma2[j])
                    del(sigma2[i])
                    b = True
                i += 1
            j += 1
        b = False
    return sigma2

# Et on fait aussi une fonction deletion_condition_theorem pour le fun

def deletionConditionTheorem(sigma,W):
    """
    Pour une liste de réflexions simples sigma représentant un mot sur S, retourne
    un mot équivalent par l'application du théorème de condition de délétion (voir
    théorème 2.7)

    Fonctionne pour les groupes de Coxeter finis (le critère de longueur théoriquement 
    vrai dans tous les groupes de Coxeter ne peut pas s'appliquer informatiquement 
    dans les groupes infinis)

    EXEMPLE ::
    
        sage: W = CoxeterGroup(['A',3])
        sage: S = W.gens()
        sage: w = [S[2],S[0],S[0],S[1],S[0],S[1],S[1],S[2],S[1],S[2]]
        sage: deletionConditionTheorem(w,W)
        [3, 2, 1, 2, 2, 3, 2, 3]
    """
    w = W.one()
    for s in sigma :
        w = w * s
    l = longueur(w,W)
    if l == len(sigma):
        return gen_to_indices(sigma,W)
    else :
        sigma2 = sigma.copy()
        for j in range(1,len(sigma2),1):
            alpha = associate_root(sigma2[j],W)
            for i in range(j):
                w = constructPartialSigma(sigma2, i, j, W)
                alphaI = w.action_on_root_indices(alpha)
                if associate_root(sigma2[i],W) == alphaI:
                    del(sigma2[j])
                    del(sigma2[i])
                    return gen_to_indices(sigma2,W)
        return -1

# Deux fonctions qui permettent de transformer un liste de générateurs en une liste de leurs indices et réciproquement

def gen_to_indices(sigma,W):
    S = W.gens()
    res = []
    for s in sigma :
        i = 0
        while s != S[i] and i < len(S) :
            i += 1
        res += [i+1]
    return res

def indices_to_gen(sigma,W):
    S = W.gens()
    res = []
    for i in sigma:
        res += [S[i-1]]
    return res

#Fonction qui étant donnée une relation de tresse, trouve toutes les occurrences de cette relation dans le mot
def indices_tresses(rel,sigma):
    res = []
    k = len(rel[0])
    for i in range(len(sigma)):
        if sigma[i:(i+k)] == rel[0] or sigma[i:(i+k)] == rel[1]:
            res += [i]
    return res

#Fonction qui étant donné une relation de tresse rel l'applique à l'indice i
def appliquer_tresse(sigma,rel,i): #Pourrait servir pour tout type de relations
    if sigma[i] == rel[0][0]:
        res = sigma[:i] + rel[1] + sigma[i+len(rel[0]):]
    else :
        res = sigma[:i] + rel[0] + sigma[i+len(rel[0]):] 
    return res

def Orbite_Tresse(sigma, W):
    """
    Pour une liste de réflexions simples sigma représentant un mot sur S, retourne
    l'ensemble des mots réduits équivalent par l'application du théorème de condition
    de délétion (voir théorème 2.7) puis de la propriété du mot (voir le théorème 3.10)

    Fonctionne pour les groupes de Coxeter finis (le critère de longueur théoriquement 
    vrai dans tous les groupes de Coxeter ne peut pas s'appliquer informatiquement 
    dans les groupes infinis)

    EXEMPLE ::
    
        sage: W = CoxeterGroup(['A',3])
        sage: S = W.gens()
        sage: w = [S[2],S[1],S[0],S[2],S[1],S[2]]
        sage: Orbite_Tresse(w,W)
        {(1, 2, 1, 3, 2, 1),
         (1, 2, 3, 1, 2, 1),
         (1, 2, 3, 2, 1, 2),
         (1, 3, 2, 1, 3, 2),
         (1, 3, 2, 3, 1, 2),
         (2, 1, 2, 3, 2, 1),
         (2, 1, 3, 2, 1, 3),
         (2, 1, 3, 2, 3, 1),
         (2, 3, 1, 2, 1, 3),
         (2, 3, 1, 2, 3, 1),
         (2, 3, 2, 1, 2, 3),
         (3, 1, 2, 1, 3, 2),
         (3, 1, 2, 3, 1, 2),
         (3, 2, 1, 2, 3, 2),
         (3, 2, 1, 3, 2, 3),
         (3, 2, 3, 1, 2, 3)}
    """
    w = gen_to_indices(reduction(sigma,W),W)
    Tmp = {tuple(w)}
    Traité = set()
    count = 0
    rels = W.braid_relations()
    for i in range(len(rels)):
        rels[i][0] = tuple(rels[i][0])
        rels[i][1] = tuple(rels[i][1])
    while len(Tmp) != 0 :
        w = Tmp.pop()
        count += 1
        print(count)
        for rel in rels :
            i = indices_tresses(rel,w)
            for j in i :
                l = appliquer_tresse(w,rel,j)
                if not(l in Tmp) and not(l in Traité):
                        Tmp.add(l)
        Traité.add(w)
    return Traité

def Orbite_Tresse_l(sigma, W):
    """
    Version peu efficace de la fonction qui calcule l'ensemble des mots réduits équivalents à sigma dans W
    Utilise des listes qui sont trop lourdes partir d'un certain rang (environ 7,5h pour sigma permutation inverse dans A5)
    """
    w = gen_to_indices(reduction(sigma,W),W)
    Tmp = [w]
    Traité = []
    count = 0
    rels = W.braid_relations()
    while len(Tmp) != 0 :
        w = Tmp.pop()
        count += 1
        print(count)
        for rel in rels :
            i = indices_tresses(rel,w)
            for j in i :
                l = appliquer_tresse(w,rel,j)
                if not(l in Tmp) and not(l in Traité):
                        Tmp.append(l)
        Traité.append(w)
    return Traité

def graphe_mots_reduits(sigma, W):
    """
    Pour une liste de réflexions simples sigma représentant un mot sur S, retourne
    l'ensemble des mots réduits équivalent par l'application du théorème de condition
    de délétion (voir théorème 2.7) puis de la propriété du mot (voir le théorème 3.10)

    Fonctionne pour les groupes de Coxeter finis (le critère de longueur théoriquement 
    vrai dans tous les groupes de Coxeter ne peut pas s'appliquer informatiquement 
    dans les groupes infinis)

    EXEMPLE ::
    
        sage: W = CoxeterGroup(['A',3])
        sage: S = W.gens()
        sage: w = [S[2],S[1],S[0],S[2],S[1],S[2]]
        sage: graphe_mots_reduits(w,W)
        Pour une liste de réflexions simples sigma représentant un mot sur S, retourne
    l'ensemble des mots réduits équivalent par l'application du théorème de condition
    de délétion (voir théorème 2.7) puis de la propriété du mot (voir le théorème 3.10)

    Fonctionne pour les groupes de Coxeter finis (le critère de longueur théoriquement 
    vrai dans tous les groupes de Coxeter ne peut pas s'appliquer informatiquement 
    dans les groupes infinis)

    EXEMPLE ::
    
        sage: W = CoxeterGroup(['A',3])
        sage: S = W.gens()
        sage: w = [S[2],S[1],S[0],S[2],S[1],S[2]]
        sage: G = graphe_mots_reduits(w,W)
        sage: G.vertices()
        [(3, 2, 1, 3, 2, 3),
         (3, 2, 3, 1, 2, 3),
         (3, 2, 1, 2, 3, 2),
         (3, 1, 2, 1, 3, 2),
         (2, 3, 2, 1, 2, 3),
         (1, 3, 2, 1, 3, 2),
         (3, 1, 2, 3, 1, 2),
         (1, 3, 2, 3, 1, 2),
         (1, 2, 3, 2, 1, 2),
         (1, 2, 3, 1, 2, 1),
         (2, 3, 1, 2, 1, 3),
         (2, 1, 3, 2, 1, 3),
         (2, 3, 1, 2, 3, 1),
         (1, 2, 1, 3, 2, 1),
         (2, 1, 3, 2, 3, 1),
         (2, 1, 2, 3, 2, 1)]
    """
    w = gen_to_indices(reduction(sigma,W),W) # on cherche d'abord un mot réduit w équivalent à sigma
    G = graphs.EmptyGraph()
    Tmp = {tuple(w)} # set qui contient les mots que l'on doit encore traiter (ie checker les relations de tresse)
    Traité = set() # set qui contient la liste des mots obtenus
    count = 0
    rels = W.braid_relations() # on récupère les relations de tresse de W
    for i in range(len(rels)): # et on les convertit en tuples pour faciliter l'opération
        rels[i][0] = tuple(rels[i][0])
        rels[i][1] = tuple(rels[i][1])
    while len(Tmp) != 0 : # quand il n'y a plus de mot à traiter (fonctionne par lemme de Matsumoto/propriété du mot)
        w = Tmp.pop() # on prend un mot au hasard
        count += 1
        print(count) # juste pour voir comment ça avance
        for rel in rels : #pour chaque relation de tresse
            i = indices_tresses(rel,w) # on calcule les emplacements de chaque occurence de la relation
            for j in i : # puis on l'applique à chaque emplacement
                l = appliquer_tresse(w,rel,j)
                G.add_edge(w,l,rel) # on ajoute l'arête correspondante au graphe dans tous les cas
                if not(l in Tmp) and not(l in Traité): # on teste l'appartenance pour voir si on a déjà rencontré ce mot ou pas
                        Tmp.add(l) # puis on l'ajoute aux mots à traiter si ce n'est pas le cas
        Traité.add(w) # on a traité w
    return G

def evacuation_path(n,T): #détermine le chemin d'évacuation de n dans le tableau T
    coord = [-1,-1]
    for i in range(len(T)):
        if T[i][-1] == n:
            coord = [i,len(T[i])-1]
    c = coord.copy()
    res = [c]
    while c[0] != 0 and c[1] != 0:
        #print(c)
        if T[c[0]][c[1]-1] > T[c[0]-1][c[1]] :
            c = [c[0],c[1]-1]
        else :
            c = [c[0]-1,c[1]]
        res.append(c)
    if c[0] == 0 :
        while T[c[0]][c[1]] != 0 and c[1] != 0:
            #print(c)
            c = [c[0],c[1]-1]
            res.append(c)
    else :
        while T[c[0]][c[1]] != 0 and c[0] != 0:
            #print(c)
            c = [c[0]-1,c[1]]
            res.append(c)
    #print(c)
    return res

def max(t): # Détermine l'élément le plus grand du tableau
    max = 0
    for i in t:
        for j in i :
            if j > max :
                max = j
    return max

def bijection_tableaux_mots(t) :
    """ Prend en argument un tableau de Young standard de forme n*(n-1)*...*1 et
    renvoie le mot équivalent au mot le plus long de An associé par la bijection de
    Edelman et Greene"""
    L = list() # transformation du tableau en listes pour effectuer les modifs
    for i in t :
        L.append(list(i))
    maxi = max(L) # détermination de l'élément max du tableau
    res = []
    while maxi != 0 : # on prend comme convention qu'une case vide contient 0, donc maxi == 0 <=> L est vide
        #print(maxi)
        path = evacuation_path(maxi,L) # on détermine son chemin d'évacuation
        #print(L)
        L[path[0][0]][path[0][1]] = 0 # on vide la case contenant le max
        for i in range(len(path)-1): # on échange chaque case pour que la case vide suive le chemin path
            L[path[i][0]][path[i][1]], L[path[i-1][0]][path[i-1][1]] = L[path[i-1][0]][path[i-1][1]], L[path[i][0]][path[i][1]]
        maxi = max(L) # on calcule le nouveau max
        res.append(path[0][1]+1) # on ajoute la permutation s_i à la liste
    res.reverse() # pas très utile pour faire directement une bijection mais elle est définie comme ça
    return tuple(res)

def max_ligne(l):
    max = l[0]
    for i in l :
        if i > max :
            max = i
    return max

def coxeter_knuth(n,T,S,c):
    if T == [] :
        T.append([n])
        S.append([c])
    else :
        i = 0
        while n != 0 :
            if i == len(T):
                T.append([n])
                S.append([c])
                #print(T)
                n = 0
            elif n >= max_ligne(T[i]) :
                T[i].append(n)
                S[i].append(c)
                #print(T)
                n = 0
            else :
                if n in T[i] and n+1 in T[i] :
                    n = n+1
                    #print(n)
                    i += 1
                else :
                    j = 0
                    while T[i][j] <= n :
                        j += 1
                    n, T[i][j] = T[i][j], n
                    #print(n)
                    i += 1

def bijection_mots_tableaux(l):
    T = []
    S = []
    for i in range(len(l)):
        coxeter_knuth(l[i],T,S,i+1)
    return StandardTableau(S)

# Tests 

def calcul_temps(w,W,s):
    t1 = time.time()
    Orbite_Tresse(w,W)
    t2 = time.time()
    s.reduced_words()
    t3 = time.time()
    print((t3 - t2))
    print((t2 - t1))
    if (t2 - t1) < (t3 - t2):
        print("super")
    else :
        print("mince")


"""
W = CoxeterGroup(['A',2])
S = W.gens()
w = S[0]*S[1]*S[0]*S[1]*S[0]
print(longueur(w,W))
print(S[0])
print(S[1])
print(w)

W = CoxeterGroup(['A',8], implementation = 'permutation')
for s in W.reflections():
    print(associate_root(s,W))

W = CoxeterGroup(['A',3])
a,b,c = W.gens()
sigma = [a,b,a,b,c,b,a,b,c,a]
print(a)
print(" ")
print(b)
print(" ")
print(c)
print(" ")
print(a*b*c)
print(" ")
print(constructPartialSigma(sigma, 1, 5, W))

# faire des tests pour reduction
print(reduction(sigma,W))
"""
"""
W = CoxeterGroup(['A',5])
S = W.gens()
w = [S[4],S[3],S[2],S[1],S[0],S[4],S[3],S[2],S[1],S[4],S[3],S[2],S[4],S[3],S[4]]
L = Orbite_Tresse(w,W)
ST = StandardTableaux([5,4,3,2,1]).list()
test = set()
count = 0
for i in ST :
    l = bijection_tableaux_mots(i)
    if l in L and l not in test:
        count += 1
        print(count)
        test.add(l)
if len(test) == len(L):
    print("gagné")
"""

"""
ST = StandardTableaux([5,4,3,2,1]).list()
for i in range(len(ST)):
    if ST[i] != bijection_mots_tableaux(bijection_tableaux_mots(ST[i])):
        print("incohérence")
    else :
        print(i+1)
"""

"""
W = CoxeterGroup(['A',5])
S = W.gens()
w = [S[4],S[3],S[2],S[1],S[0],S[4],S[3],S[2],S[1],S[4],S[3],S[2],S[4],S[3],S[4]]
L = Orbite_Tresse(w,W)
for i in range(len(L)):
    l = L.pop()
    if l != bijection_tableaux_mots(bijection_mots_tableaux(l)):
        print("incohérence")
    else :
        print(i+1)
"""