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

def associate_root(sigma,W): # Pour sigma reflexion, renvoie la racine positive qui lui est associée
    for alpha in range(len(W.positive_roots())):
        if W.roots()[sigma.action_on_root_indices(alpha)] == - W.positive_roots()[alpha] :
            return W.positive_roots()[alpha]
    return -1

# On écrit une fonction permettant la construction des
# applications partielles sigma[i+1] ... sigma[j-1]
def constructPartialSigma(sigma, i, j, W):
    w = W.one()
    for s in sigma[i+1:j]:
        w = w * s
    return w

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