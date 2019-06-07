import numpy as np

def match(a, b):
  b_set = set(b)
  b_match = [i for i, v in enumerate(a) if v in b_set]
  a_set = set(a)
  a_match = [i for i, v in enumerate(b) if v in a_set]
  
  a_match = np.asarray(a_match)
  b_match = np.asarray(b_match)
  
  a_match2=np.argsort(a[b_match])
  b_match2=np.argsort(b[a_match])
  
  return b_match[a_match2],a_match[b_match2]