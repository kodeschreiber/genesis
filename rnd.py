import math
import hashlib
import typing
import yaml
import multiprocessing as mp
from dataclasses import dataclass

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

  def iter(self, **args):
    for i in range(self.range(**args)):
      yield i

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
rstar = (10, 1000)
earth_vol = seedobj.frange(0.55, 0.90)
water_vol = seedobj.frange(0.3, 0.35)
air_vol = seedobj.frange(0.01, 0.90)
fire_vol = seedobj.frange(0.01, 0.30)
stars = list()  # (float: intensity, int: radius, )
star_odds = [ 0.85, 0.25, 0.05 ]
steps = 100
precipes = None
lrecipes = None


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
    self.x = int(self.__pos / self.__length)
    self.y = self.__pos % self.__length
    if not self.__fully_init:
      res = (self.x, self.y)
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
    self.lycan = 0.0
    self.weeds = 0.0
    self.dirt = 0.0
    self.grass = 0.0
    self.mud = 0.0
    self.cloud = 0.0
    self.clear = 0.0
    self.clearwater = 0.0
    self.murky = 0.0
    self.alginated = 0.0
    self.magmatic = 0.0
    self.scorched = 0.0
    self.snow = 0.0
    self.ice = 0.0
    self.fog = 0.0
    self.smoke = 0.0
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

@dataclass
class Celestial:
  cx: int
  cy: int
  lumosity: float
  peak: float
  speed: float
  arc: float
  offset: float
  direction: bool

  def move(self, dx):
    dire = 1 if self.direction else -1
    self.cx = dx * self.speed
    self.cy = dire * (arc * dx) + offset

  def rads_at(pt):
    dist = math.sqrt((pt[0] - ctr[0])**2 + (pt[1] - ctr[1])**2)
    return ((-0.1 * dist**2 + self.peak) / self.peak) * self.lumosity

  def sr_star(seedobj, width, length):
    return Celestial( seedobj.range(width),
                      seedobj.range(length),
                      seedobj.frange(0.05, 0.95),
                      seedobj.frange(0.10, 0.80) * int((width + length)/2),
                      seedobj.frange(0.001, 1.0) * 10,
                      seedobj.frange(0.001, 1.0) * 10,
                      seedobj.range(width/2),
                      seedobj.odds(0.5),
                    )

@dataclass
class Recipe:
  inputs: typing.Dict
  outputs: typing.Dict

  def volume_check(self):
    return round(sum(self.inputs.values()), 5) == round(sum(self.outputs.values()), 5)

  def __can_perform__(self, pt):
    for k,v in self.inputs.items():
      if getattr(pr, k) < v:
        return False
    return True

  def attempt(self, pt):
    if self.__can_perform__(pt):
      for k,v in self.inputs.items():
        setattr(pt, k, getattr(pt, k) - v)
      for k,v in self.outputs.items():
        setattr(pt, k, getattr(pt, k) + v)

  def dryrun(self, pt):
    pass

  def total(self):
    return round(sum(self.inputs.values()), 5) * 2


@dataclass
class LateralRecipe:
  origin: Recipe
  destination: Recipe
  radius: float

  def move(self):
    pass


def give_me_a_sine(seedobj, frequency, sample_rate, variance):
  return [ Sine.randinit(seedobj, frequency, sample_rate) for i in range(variance) ]

def you_are_my_sumsine(sinelist, pt):
  return sum([ j.calc(pt) for j in sinelist ])

def level_with_me(grid, attr, top):
  Z = [ getattr(pt, attr) for pt in grid ]
  low = min(Z)
  tmax = top/float(max(Z) + low)
  for pt in grid:
    setattr(pt, attr, int((z + low) * tmax))

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

with open('recipes.yaml', 'r') as yfile:
  data = yaml.load(yfile, yaml.Loader)
  precipes = [ Recipe(**rec) for rec in data['point'] ]
  lrecipes = [ LateralRecipe(**rec) for rec in data['lateral'] ]

# At first, there was a plane, as flat to the horizons
grid = IterGrid(width, length)
for pt in grid:
  grid.set(pt, Point())

# Then, by the shift of the earth, great mounds and valleys made their appearance
X = give_me_a_sine(seedobj, frequency, sample_rate, variance)
Y = give_me_a_sine(seedobj, frequency, sample_rate, variance)
for pt in grid:
  pt.smooth = you_are_my_sumsine(grid.x, X) * you_are_my_sumsine(grid.y, Y)
level_with_me(grid, 'smooth', earth_vol)

# From the ground was water ammased at all the lowest points
for pt in grid:
  if pt.earth() < water_vol:
    pt.clearwater = water_vol - pt.earth()

# As the surface settles, extra-planar gasses decend upon the region
for pt in grid
  if pt.earth() + pt.water() < air_vol:
    pt.clear = air_vol - (pt.earth() + pt.water())

# From below, magmatic forces are driven forth
X = give_me_a_sine(seedobj, frequency, sample_rate, variance)
Y = give_me_a_sine(seedobj, frequency, sample_rate, variance)
for pt in grid:
  pt.magmatic = you_are_my_sumsine(grid.x, X) * you_are_my_sumsine(grid.y, Y)
level_with_me(grid, 'magmatic', fire_vol)

# From the heavens were the celestials strewn, casting their light over all
for val in star_odds:
  if seedobj.odds(val):
    stars.append(Celestial.sr_star(seedobj, width, length))

# From the heavens the Alchemist conjured massive hydroplanes to blanket the sky. From them
# would many droplets of water precipitate upon the land, filing their way to it's lowest points
settled = False
while not settled:
  settled = True
