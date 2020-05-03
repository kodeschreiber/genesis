import math
import hashlib
import typing
import multiprocessing as mp

class Seeder:
  _min = 2
  _max = 8
  _hash_length = 64
  def __init__(self, seed):
    self._ptr = 0
    self._val = seed
    self._hash_()
    self._count = (self._val[0] % (Seeder._max-Seeder._min)) + Seeder._min

  def _hash_(self):
    self._val = hashlib.sha256(self._val).hexdigest()

  def _get_value_(self):
    if (self._ptr + self._count) > Seeder._hash_length:
      self._hash_()
      self._ptr = 0
    return self._val[self._ptr:self._ptr+self._count]

  def float(self):
    return int(self._get_value_(), 16) / float(2**(self._count * 4))

  def range(self, **args):
    min = 0
    max = 0
    if len(args) == 1:
      max = args[0]
    elif len(args) == 2:
      min,max = args
    else:
      max = 2**(self._count * 4)
    return (int(self._get_value_(), 16) % (max-min)) + min

  def frange(self, **args):
    prec = -1
    for x in args:
      cur = len(str(x).split('.')[1])
      if cur > prec:
        prec = cur
    return self.range([ int(x*(10**prec)) for x in args ])/float(10**prec)

  def polar(self):
    return -1 if (self._get_value_ % 2) == 0 else 1

  def fpolar(self):
    return self.polar() * self.float()

  def ipolar(self, **args):
    return self.polar() * self.range(**args)

  def odds(self, float_success):
    return self.float() <= float_success

# The Alchemist's Genesis

# As with many a histories of geometrical nature, it began with a compass and a square
width = 1000
length = 1000

# Each literal had measure
seedobj = Seeder("Abacab")
scale = 10
sample_rate = 20.0
frequency = 1
variance = 30
flatness = 0.0001
water_vol = seedobj.frange(0.3, 0.35)
air_vol = seedobj.frange(0.01, 0.90)
mantle_vol = seedobj.frange(0.01, 0.30)

# And each modus operandi was enuciated in runic fashion within the lexicon of the Alchemist;
# for change is the nature of all, some of the greatest, some of the slight

class IterGrid:
  def __init__(self, width, length):
    self.__pos = 0
    self.__fully_init = False
    self.__width = width
    self.__length = length
    self.__grid_data = list()
    self.__grid_len = width * length
    for i in range(width*length):
      self.__grid_data.append(None)

  def __iter__(self):
    self.__pos = 0
    return self

  def __next__(self):
    res = None
    if not self.__fully_init:
      res = (int(self.__pos / self.__length), self.__pos % self.__length)
    else:
      res = self.__grid_data[self.__pos]

    if self.__pos < self.__grid_len:
      self.__pos += 1
    else:
      if not self.__fully_init:
        self.__fully_init = True
      raise StopIteration
    return res

  def set(self, pt, val):
    print(pt, (pt[0]*self.__length)+pt[1])
    self.__grid_data[(pt[0]*self.__length)+pt[1]] = val

class Point:
  capacity = 256
  def __init__(self):
    self.igneous = 0.0
    self.smooth = 0.0
    self.boulder = 0.0
    self.cobbled = 0.0
    self.gravel = 0.0
    self.sand = 0.0
    self.lycanned = 0.0
    self.weeded = 0.0
    self.dirt = 0.0
    self.grassy = 0.0
    self.mud = 0.0
    self.cloud = 0.0
    self.clear = 0.0
    self.clearwater = 0.0
    self.murky = 0.0
    self.alginated = 0.0
    self.magmatic = 0.0
    self.scorched = 0.0
    self.snow = 0.0
    self.iced = 0.0
    self.fogged = 0.0
    self.smogged = 0.0
    self.rain = 0.0
    self.wind = 0.0
    self.hail = 0.0
    self.solar = 0.0
    self.lava = 0.0
    self.blaze = 0.0

  def earth(self):
    return self.igneous + self.smooth + self.boulder + self.cobbled + self.gravel + self.sand + self.lycanned + self.weeded + self.dirt + self.grass + self.mud + self.lava

  def water(self):
    return self.alginated + self.murky + self.clearwater

  def air(self):
    return self.clear + self.fogged + self.cloud + self.wind + self.hail

  def fire(self):
    return self.magmatic + self.scorched + self.solar = self.blaze

  def total(self):
    return self.earth() + self.water() + self.air() + self.fire()

class Sine:
  def __init__(self, scalar, frequency, sample_rate, offset):
    self.scalar = scalar
    self.frequency = frequency
    self.sample_rate = sample_rate
    self.offset = offset

  def calc(self, x):
    return self.scalar * math.sin(2*math.pi * self.frequency * (x/self.sample_rate)) + self.offset

  def calc_range(self, min, max):
    for i in range(min, max):
      yield self.calc(i)

  @staticmethod
  def randinit(seedobj, frequency, sample_rate):
    return Sine(seedobj.fpolar(), seedobj.ripolar(frequency), sample_rate, seedobj.fpolar())

def calc_mesh(seedobj, variance, frequency, sample_rate, width, length):
  ret = list
  X = [ Sine.randinit(seedobj, frequency, sample_rate) for i in range(variance) ]
  Y = [ Sine.randinit(seedobj, frequency, sample_rate) for i in range(variance) ]
  for xpt in range(width):
    for ypt in range(length):
      ret.append((xpt, ypt, sum([ j.calc(xpt) for j in X ]) * sum([ k.calc(ypt) for k in Y ])))

def level(zmesh, top):
  Z = [ i[2] for i in zmesh ]
  low = min(Z)
  tmax = top/float(max(Z) + low)
  for x, y, z in zmesh:
    yield (x, y, int((z + low) * tmax))

def meshran(xmax, ymax):
  for i in range(xmax):
    for j in range(ymax):
      yeild (i, j)

def get_stars(seedobj):
  pass

def prox(pt, rad, dim):
  pts = list()
  spt = (pt[0] - rad, pt[1] - rad)
  ept = (pt[0] - rad, pt[1] - rad)
  if spt[0] < 0:
    spt[0] = 0
  if spt[1] < 0:
    spt[1] = 0
  if ept[0] > dim[0]:
    ept[0] = dim[0]-1
  if ept[1] > dim[1]:
    ept[1] = dim[1]-1
  for x in range(ept[0]-spt[0]):
    for y in range(ept[1]-spt[1]):
      if math.sqrt((x - pt[0])**2 + (y - pt[0])**2) <= rad:
        pts.append((x,y))
  return pts

def fshift(pts):


# At first, there was a plane, as flat to the horizons
grid = list()
for i in range(width):
  grid[i] = list()
  for j in range(length):
    grid[i][j] = Point()

# Then, by the shift of the earth, great mounds and valleys made their appearance
pts = calc_mesh(seedobj, variance, frequency, sample_rate, width, length)
for x, y, z in level(pts, 1):
  grid[x][y].smooth = z

# From the ground was water ammased at all the lowest points
for x, y in meshran(width, length):
  if grid[x][y].earth() < water_vol:
    grid[x][y].clearwater = water_vol - grid[x][y].earth()

# As the surface settles, extra-planar gasses decend upon the region
for x, y in meshran(width, length):
  if grid[x][y].earth() + grid[x][y].water() < air_vol:
    grid[x][y].clear = air_vol - (grid[x][y].earth() + grid[x][y].water())

# From below, magmatic forces are driven forth
pts = calc_mesh(seedobj, variance, frequency, sample_rate, width, length)
for x, y, z in level(pts, mantle_vol):
  grid[x][y].magmatic = z

# From the heavens the Alchemist conjured massive hydroplanes to blanket the sky. From them
# would many droplets of water precipitate upon the land, filing their way to it's lowest points
settled = False
while not settled:
  settled = True
