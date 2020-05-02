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

  def range(**args):
    min = 0
    max = 0
    if len(args) == 1:
      max = args[0]
    elif len(args) == 2:
      min,max = args
    else:
      max = 2**(self._count * 4)
    return (int(self._get_value_(), 16) % (max-min)) + min

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
water_vol = seedobj.range(3, 35)/100.0
air_vol = seedobj.range(1, 90)/100.0

# And each modus operandi was enuciated in runic fashion within the lexicon of the Alchemist;
# for change is the nature of all, some of the greatest, some of the slight

class Point:
  capacity = 256
  def __init__(self):
    igneous = 0.0
    smooth = 0.0
    boulder = 0.0
    cobbled = 0.0
    gravel = 0.0
    sand = 0.0
    lycanned = 0.0
    weeded = 0.0
    dirt = 0.0
    grassy = 0.0
    mud = 0.0
    cloud = 0.0
    clear = 0.0
    clearwater = 0.0
    murky = 0.0
    alginated = 0.0
    magmatic = 0.0
    scorched = 0.0
    snow = 0.0
    iced = 0.0
    fogged = 0.0
    smogged = 0.0
    rain = 0.0
    wind = 0.0
    hail = 0.0
    solar = 0.0

  def earth(self):
    return self.igneous + self.smooth + self.boulder + self.cobbled + self.gravel + self.sand + self.lycanned + self.weeded + self.dirt + self.grass + self.mud

  def water(self):
    return self.alginated + self.murky + self.clearwater

  def air(self):
    return self.clear + self.fogged + self.cloud + self.wind + self.hail

  def fire(self):
    return self.magmatic + self.scorched + self.solar

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
    if grid[x][y].total() > 1:
      grid[x][y].clear = 1 - (grid[x][y].earth() + grid[x][y].water())

# From the heavens the Alchemist conjured massive hydroplanes to blanket the sky. From them
# would many droplets of water precipitate upon the land, filing their way to it's lowest points
settled = False
while not settled:
  settled = True
