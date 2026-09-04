R = SR   

sigma = diagonal_matrix(R, [
    1,1,1,-1,-1,-1,
     1, 1, 1, -1, -1, -1,
    1,
     1
])


iota = -diagonal_matrix(R, [
    -1,-1,-1,-1,-1,-1,
     1, 1, 1, 1, 1, 1,
    -1,
     1
])

Z6 = zero_matrix(R,6)
I6 = identity_matrix(R,6)

rho = block_matrix(R, [
    [ Z6, -I6,              zero_matrix(R,6,1), zero_matrix(R,6,1)],
    [ I6,  Z6,              zero_matrix(R,6,1), zero_matrix(R,6,1)],
    [ zero_matrix(R,1,6), zero_matrix(R,1,6), matrix(R,[[0]]), matrix(R,[[-1]])],
    [ zero_matrix(R,1,6), zero_matrix(R,1,6), matrix(R,[[1]]), matrix(R,[[0]])]
])

g_diag = exp(matrix(R, [
    [ pi*I/4, 0],
    [0,-pi*I/4]
]))

g_off = exp(matrix(R, [
    [0,      pi*I/4],
    [pi*I/4, 0]
]))

g_offskew = matrix([[1/sqrt(2),1/sqrt(2)],[-1/sqrt(2),1/sqrt(2)]])


def identity14():
    return identity_matrix(R,14)


idx = {
"xi1":0,"xi2":1,"xi3":2,
"c1":3,"c2":4,"c3":5,
"eta1":6,"eta2":7,"eta3":8,
"d1":9,"d2":10,"d3":11,
"xi":12,"eta":13
}



def compute(g):
    g_inv = g.inverse()
    g_Cartan = g_inv.transpose()

    def apply_pair(M,i,j):
        a,b,c,d = g.list()
        M[i,i]=a
        M[i,j]=b
        M[j,i]=c
        M[j,j]=d

    def apply_pair2(M,i,j):
        a,b,c,d = g_Cartan.list()
        M[i,i]=a
        M[i,j]=b
        M[j,i]=c
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

    beta = A1*A2*A3

    if g == g_diag:
        beta *= sqrt(2) * (1 + I) / 2
    elif g == g_off:
        beta *= 2*sqrt(2)
    elif g == g_offskew:
        beta *= 2*sqrt(2)

    cocycle = beta.inverse()*beta.conjugate()

    print("beta =")
    print(beta)

    print("cocycle =")
    print(beta.inverse()*beta.conjugate())

# beta.conjugate()*beta.inverse() も考えてみようとしたけどこの場合一致しました。

    print("beta*sigma*beta^{-1}=")
    print(beta*sigma*beta.inverse())

    print("beta*rho*beta^{-1}=")
    print(beta*rho*beta.inverse())

    print("beta*iota*beta^{-1}=")
    print(beta*iota*beta.inverse())

    print("beta*rho*iota*beta^{-1}=")
    print(beta*rho*iota*beta.inverse())



    if g == g_offskew:
        print("cocycle = 1")
        print(cocycle == identity14())
        print("beta*sigma*beta^{-1}=sigma")
        print(beta*sigma*beta.inverse()==sigma)
        print("beta*rho*beta^{-1}=rho")
        print(beta*rho*beta.inverse()==rho)
        print("beta*iota*beta^{-1}=iota*rho")
        print(beta*iota*beta.inverse()== iota*rho)
        print("beta*rho*iota*beta^{-1}=iota")
        print(beta*rho*iota*beta.inverse()==iota)


    elif g == g_off:
        print("cocycle = -sqrt{-1}*iota*rho")
        print(cocycle == -I*iota*rho)
        print("beta*sigma*beta^{-1}=sigma")
        print(beta*sigma*beta.inverse()==sigma)
        print("beta*rho*beta^{-1}=-sqrt{-1}*iota")
        print(beta*rho*beta.inverse()==-I*iota)
        print("beta*iota*beta^{-1}=-sqrt{-1}*rho")
        print(beta*iota*beta.inverse()==-I*rho)
        print("beta*rho*iota*beta^{-1}=rho*iota")
        print(beta*rho*iota*beta.inverse()==rho*iota)

    elif g == g_diag:
        print("cocycle= iota")
        print(cocycle == iota)
        print("beta*sigma*beta^{-1}=sigma")
        print(beta*sigma*beta.inverse()==sigma)
        print("beta*rho*beta^{-1}=-sqrt{-1}*iota*rho")
        print(beta*rho*beta.inverse()==-I*iota*rho)
        print("beta*iota*beta^{-1}=iota")
        print(beta*iota*beta.inverse()==iota)
        print("beta*rho*iota*beta^{-1}=sqrt{-1}*rho")
        print(beta*rho*iota*beta.inverse()==I*rho)

    return beta


print("=== beta_1 ===")
beta_offskew = compute(g_offskew)

print("=== beta_2 ===")
beta_off = compute(g_off)


print("=== beta_3 ===")
beta_diag = compute(g_diag)

