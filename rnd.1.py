import math
import hashlib
import typing

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
water_vol = seedobj.range(3, 35)

# And each modus operandi was enuciated in runic fashion within the lexicon of the Alchemist;
# for change is the nature of all, some of the greatest, some of the slight

class Point:
  capacity = 256
  earth_delta = 12
  earth_types = {
    'igneous': (0.15, 'igneous'),  # Type: (odds of moving, replaced with); inverse used for water moving odds
    'smooth': (0.01, 'boulder'),
    'boulder': (0.10, 'boulder'),
    'cobbled': (0.21, 'cobbled'),
    'gravel': (0.43, 'gravel'),
    'sand': (0.78, 'sand'),
    'lycanned': (0.01, 'boulder'),
    'weeded': (0.21, 'cobbled'),
    'dirt': (0.65, 'mud'),
    'grassy': (0.18, 'mud'),
    'mud': (0.82, 'mud')
  }
  water_types = [
    'cloud',
    'clear',
    'murky',
    'alginated'
  ]
  air_types = [
    'clear',
    'windy',
    'foggy'
  ]
  fire_types = [
    'sunny',
    'magmatic',
    'blazed'
  ]

  def __init__(self):
    self.earth = 0
    self.earth_type = None
    self.water = 0
    self.water_type = None
    self.air = 0
    self.air_type = None
    self.fire = None
    self.fire_type = None
    self._water_passed = 0

  def humidify(self, amount=1):
    self.water += 1
    self._water_passed += 1


  def drain(self, amount=1):

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

# At first, there was a plane, as flat to the horizons
grid = list()
for i in range(width):
  grid[i] = list()
  for j in range(length):
    grid[i][j] = Point()

# Then, by the shift of the earth, great mounds and valleys made their appearance
pts = calc_mesh(seedobj, variance, frequency, sample_rate, width, length)
for x, y, z in level(pts, Point.capacity):
  grid[x][y].earth = z
  grid[x][y].water = water_vol

# From the heavens the Alchemist conjured massive hydroplanes to blanket the sky. From them
# would many droplets of water precipitate upon the land, filing their way to it's lowest points
settled = False
while not settled:
  settled = True
