#!/usr/bin/env sage

# author : Lucas POUILLART

# On va commencer par redéfinir la fonction de réduction pour les groupes de Weyl, afin de l'appliquer
# aux groupes de Coxeter en général :

# On utilise le fait que dans sage, la fonction roots() renvoie une liste contenant d'abord les racines
# positives puis les négatives
def longueur(sigma, W) : # à tester plus
    pi = W.positive_roots()
    minus_pi = W.roots()[len(pi):]
    res = [] #liste qui contiendra l'intersection de pi et de sigma ^-1 (-pi)
    #boucle qui teste pour tout alpha dans pi si sigma(alpha) dans minus_pi, et si oui l'ajoute à res
    for alpha in range(len(pi)) :
        beta = W.roots()[sigma.action_on_root_indices(alpha)]
        if beta in minus_pi :
            res += [alpha]
    return len(res)

def associate_root(sigma,W): # Pour sigma reflexion, renvoie l'indice de la racine positive qui lui est associée
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
    res = sigma.copy()
    if res[i:i+len(rel[0])] == rel[0]:
        res[i:i+len(rel[0])] = rel[1]
    else :
        res[i:i+len(rel[1])] = rel[0]
    return res

def Orbite_Tresse(sigma, W):
    w = gen_to_indices(reduction(sigma,W),W)
    Tmp = [w]
    Traité = []
    rels = W.braid_relations()
    while Tmp != [] :
        w = Tmp[0]
        print(w)
        for rel in rels :
            #print(rel)
            i = indices_tresses(rel,w)
            for j in i :
                #print(j)
                l = appliquer_tresse(w,rel,j)
                if l not in Tmp and l not in Traité:
                    Tmp += [l]
        Traité += [Tmp[0]]
        del(Tmp[0])
    return Traité

# Tests 
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