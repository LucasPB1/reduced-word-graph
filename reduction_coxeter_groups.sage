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

# Tests 
W = CoxeterGroup(['A',2])
S = W.gens()
w = S[0]*S[1]*S[0]*S[1]*S[0]
print(longueur(w,W))
print(S[0])
print(S[1])
print(w)