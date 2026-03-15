R = SR   # symbolic ring (complex numbers OK)


g = exp(matrix([[pi*I/4,0],[0,-pi*I/4]]))


g_inv = g.inverse()

def identity14():
    return identity_matrix(R,14)


idx = {
"xi1":0,"xi2":1,"xi3":2,
"c1":3,"c2":4,"c3":5,
"eta1":6,"eta2":7,"eta3":8,
"d1":9,"d2":10,"d3":11,
"xi":12,"eta":13
}


def apply_pair(M,i,j):
    a,b,c,d = g.list()
    M[i,i]=a
    M[i,j]=b
    M[j,i]=c
    M[j,j]=d


def apply_pair2(M,i,j):
    a,b,c,d = g_inv.list()
    M[i,i]=a
    M[i,j]=c
    M[j,i]=b
    M[j,j]=d


A1 = identity14()
apply_pair(A1, idx["xi1"], idx["eta"])
apply_pair(A1, idx["xi"], idx["eta1"])
apply_pair(A1, idx["eta2"], idx["xi3"])
apply_pair(A1, idx["eta3"], idx["xi2"])
apply_pair2(A1, idx["c1"], idx["d1"])


A2 = identity14()
apply_pair(A2, idx["xi2"], idx["eta"])
apply_pair(A2, idx["xi"], idx["eta2"])
apply_pair(A2, idx["eta1"], idx["xi3"])
apply_pair(A2, idx["eta3"], idx["xi1"])
apply_pair2(A2, idx["c2"], idx["d2"])


A3 = identity14()
apply_pair(A3, idx["xi3"], idx["eta"])
apply_pair(A3, idx["xi"], idx["eta3"])
apply_pair(A3, idx["eta1"], idx["xi2"])
apply_pair(A3, idx["eta2"], idx["xi1"])
apply_pair2(A3, idx["c3"], idx["d3"])


delta0 = A1*A2*A3
delta = sqrt(2)*(1 + I)/2*delta


L = matrix(R,14,14)

L[idx["xi1"],idx["eta1"]] = 1
L[idx["xi2"],idx["eta2"]] = 1
L[idx["xi3"],idx["eta3"]] = 1
L[idx["c1"],idx["d1"]] = 1
L[idx["c2"],idx["d2"]] = 1
L[idx["c3"],idx["d3"]] = 1

L[idx["eta1"],idx["xi1"]] = -1
L[idx["eta2"],idx["xi2"]] = -1
L[idx["eta3"],idx["xi3"]] = -1
L[idx["d1"],idx["c1"]] = -1
L[idx["d2"],idx["c2"]] = -1
L[idx["d3"],idx["c3"]] = -1

L[idx["xi"],idx["eta"]] = 1
L[idx["eta"],idx["xi"]] = -1



result = delta.inverse() * L * delta

print("alpha1 =")
print(A1)
print("alpha2 =")
print(A2)
print("alpha3 =")
print(A3)

cocycle = delta.inverse() * delta.conjugate()

print("delta:=(1+sqrt{-1})/sqrt(2)alpha1 * alpha2 * alpha3 =")

print(delta)

print("delta^{-1} * lambda * delta =")

print(result)

print("cocycle =")

print(cocycle)