import copy

def min_S_value(D, u):
    min_S, min_i, min_j = float("inf"),-1,-1
    for k in D:
        for l in D[k]:
            if l!=k:
                crit = D[k][l] - u[k] - u[l]
                if crit < min_S:
                    min_S = crit
                    min_i = k
                    min_j = l
    return (min_i, min_j)

def neighbor_join(D):
    D = copy.deepcopy(D)

    T = {k:{} for k in D}
    r = len(D)
    while len(D)>1:
      # Find clusters minimizing D(ci, cj) - u(ci) - u(cj)
      u = {k : sum(D[k].values()) / (len(D) - 2) for k in D} if len(D) > 2 else {k : sum(D[k].values()) for k in D}
      i, j = min_S_value(D, u)

      #Compute D(r, m) for every other cluster m and add to D
      D[r] = {}
      for m in D:
        if m not in (i, j, r):
          D[r][m] = D[m][r] = 0.5*(D[i][m] + D[j][m] - D[i][j])

      # Add new vertex r to T with edges to ci and cj
      T[r] = {}
      T[i][r] = T[r][i] = 0.5*(D[i][j] + u[i] - u[j])
      T[j][r] = T[r][j] = 0.5*(D[i][j] + u[j] - u[i])

      # remove ci and cj from D
      del D[i]
      del D[j]

      for m in D:
        if i in D[m]:
          del D[m][i]
        if j in D[m]:
          del D[m][j]

      r += 1

    return T


def relaxed_neighbor_join(D):
    D = copy.deepcopy(D)

    T = {k:{} for k in D}
    r = len(D)
    while len(D)>2:
      n = len(D)
      i, j, ui, uj = None, None, None, None
      for a in D:
        S_a = {D[a][k] - (sum(D[a].values()) + sum(D[k].values()) - 2*D[a][k]) / (n - 2) for k in D if k != a}
        for b in D:
          if a == b:
            continue
          S_b = {D[b][k] - (sum(D[b].values()) + sum(D[k].values()) - 2*D[b][k]) / (n - 2) for k in D if k != b}


          ui = (sum(D[a].values()) - D[a][b]) / (n - 2)
          uj = (sum(D[b].values()) - D[a][b]) / (n - 2)
          S_ab = D[a][b] - ui - uj

          if S_ab <= min(S_a) and S_ab <= min(S_b):
            i, j = a, b
            break

        if i != None and j != None:
          break
            
      T[r] = {}
      T[i][r] = T[r][i] = 0.5*(D[i][j] + ui - uj)
      T[j][r] = T[r][j] = 0.5*(D[i][j] + uj - ui)
      
      #Compute D(r, m) for every other cluster m and add to D
      D[r] = {}
      for m in D:
        if m not in (i, j, r):
          D[r][m] = D[m][r] = 0.5*(D[i][m] + D[j][m] - D[i][j])

      # remove ci and cj from D
      del D[i]
      del D[j]

      for m in D:
        if i in D[m]:
          del D[m][i]
        if j in D[m]:
          del D[m][j]

      r += 1

    i, j = D.keys()
    T[r] = {}
    T[i][r] = T[r][i] = T[j][r] = T[r][j] = 0.5*D[i][j]
    
    return T