from copy import deepcopy
import numpy as np

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
    D = deepcopy(D)

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
    D = deepcopy(D)

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


def rapid_neighbor_join(D, idx2node):
    # initialize data structures
    D = deepcopy(D)
    idx2node = deepcopy(idx2node)
    I = {}
    S = {}
    for row in D:
      I[row] = np.argsort(D[row])
      S[row] = np.take_along_axis(D[row], I[row], axis=0)

    merged_nodes = set()
    T = {k: {} for k in idx2node}

    m_idx = len(D)
    node_count = len(D)

    while node_count > 2:
      u = {}
      for row in D:
        if row in merged_nodes:
          continue

        u[row] = 0
        for col in range(len(D)):
          if idx2node[col] not in merged_nodes:
            u[row] += D[row][col] / (node_count - 2)

      umax = max(u.values())
      qmin, i, j = float("inf"), -1, -1

      for r in range(len(S)):
        r_node = idx2node[r]
        if r_node in merged_nodes:
          continue

        for c in range(len(S[r_node])):
          curr_j = I[r_node][c]
          j_node = idx2node[curr_j]
          if j_node in merged_nodes or r == curr_j:
            continue
          
          curr_S = S[r_node][c]
          if curr_S - u[r_node] - umax > qmin:
            break
          
          curr_Q = D[r_node][curr_j] - u[r_node] - u[j_node]
          if curr_Q < qmin:
            qmin, i, j = curr_Q, r, curr_j

      i_node = idx2node[i]
      j_node = idx2node[j]

      m = str(m_idx)
      T[m] = {}
      T[i_node][m] = T[m][i_node] = 0.5*(D[i_node][j] + u[i_node] - u[j_node])
      T[j_node][m] = T[m][j_node] = 0.5*(D[i_node][j] + u[j_node] - u[i_node])


      merged_nodes.update([i_node, j_node])
      idx2node = np.append(idx2node, m)
      D[m] = np.zeros(len(idx2node))

      for r in range(len(D)):
        r_node = idx2node[r]
        if r_node not in merged_nodes and r_node != m:                  
          D_xr = 0.5*(D[i_node][r] + D[j_node][r] - D[i_node][j])     # distance value for new node to other nodes
          D[r_node] = np.append(D[r_node], D_xr)                      # add new distance to end of each row
          D[m][r] = D_xr                                    # update new row for D with new distance value
            
          insert_idx = np.searchsorted(S[r_node], D_xr)          # index for where to put new node distance in existing array
          I[r_node] = np.insert(I[r_node], insert_idx, m_idx)             # insert index for new node into row for I
          S[r_node] = np.insert(S[r_node], insert_idx, D_xr)          # insert new node distance into corresponding row for new S array

      I[m] = np.argsort(D[m])
      S[m] = np.take_along_axis(D[m], I[m], axis=0)

      m_idx += 1
      node_count -= 1

    i, j = [k for k in D if k not in merged_nodes]
    m = str(m_idx)
    T[m] = {}
    T[i][m] = T[m][i] = T[j][m] = T[m][j] = 0.5*D[i][int(j)]

    return T