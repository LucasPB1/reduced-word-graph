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

# Tests 

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