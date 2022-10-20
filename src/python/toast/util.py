def Linesearch(x0, d, s0, p0, func):
  sl = 0
  pl = p0
  sh = s0
  x = x0 + d * sh
  ph = func(x)

  if ph < pl:
    sm = sh
    pm = ph
    sh = sh * 2
    x = x0 + d * sh
    ph = func(x)
    while ph < pm:
      sl = sm
      pl = pm
      sm = sh
      pm = ph
      sh = sh * 2
      x = x0 + d * sh
      ph = func(x)

  else:
    sm = sh / 2
    x = x0 + d * sm
    pm = func(x)
    while pm > pl:
      sh = sm
      ph = pm
      sm = sm / 2
      x = x0 + d * sm
      pm = func(x)

  if ph < pm:
    pmin = pm
    s = sm
  else:
    a = ((pl - ph) / (sl - sh) - (pl - pm) / (sl - sm)) / (sh - sm)
    b = (pl - ph) / (sl - sh) - a * (sl + sh)
    s = -b / (2 * a)
    x = x0 + d * s
    pmin = func(x)
    if pmin > pm:
      s = sm
      pmin = pm

  return (s, pmin)