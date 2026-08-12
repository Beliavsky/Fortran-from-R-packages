

# sk stores the existing sequence of subset.
def btsubsetSum(v, lb, ub, uniqueOut = True, float2int64 = True):
  if type(lb) is not list: lb = [lb]; ub = [ub]
  if float2int64: v = [int(x) for x in v]
  initialLB = [lb[i] for i in range(len(lb))]
  csum = [0] * (len(lb) + 1)
  for i in range(1, len(csum) - 1): csum[i] = csum[i - 1] + v[lb[i - 1]]
  rst = set()
  if not uniqueOut: rst = []
  # atom = 0
  while True:
    if lb[-1] >= ub[-1]: # Rearrange the pointers.
      csum[-1] = csum[-2] + v[lb[-1]]
      if uniqueOut: rst.add(csum[-1])
      else: rst.append(csum[-1])
      # Find the first pointer that has not reached the upper bound:
      i = len(ub) - 1
      while i >= 0 and lb[i] >= ub[i]: i -= 1
      if i < 0: break
      lb[i] += 1
      csum[i + 1] += v[lb[i]] - v[lb[i] - 1]
      i += 1
      while i < len(ub) - 1:
        lb[i] = max(initialLB[i], lb[i - 1] + 1)
        csum[i + 1] = csum[i] + v[lb[i]]
        i += 1
      lb[i] = max(initialLB[i], lb[i - 1] + 1)
    else:
      csum[-1] = csum[-2] + v[lb[-1]]
      if uniqueOut: rst.add(csum[-1])
      else: rst.append(csum[-1])
      lb[-1] += 1
    # atom += 1
    # print("atom = ", atom, "\n")
  return list(rst)




def NofSums(lb, ub):
  if type(lb) is not list: lb = [lb]; ub = [ub]
  if sum(map(lambda x: abs(x[0] - x[1]), zip(lb, ub))) == 0: return 1;
  current = list(range(1, ub[0] - lb[0] + 1 + 1))
  for i in range(1, len(lb)):
    m = min(lb[i] - lb[i - 1], len(current))
    thenext = [current[m - 1]]
    for j in range(lb[i] + 1, ub[i] + 1):
      m = min(j - lb[i - 1], len(current))
      thenext.append(current[m - 1] + thenext[-1])
    thenext, current = current, thenext
  return current[-1]
      

































