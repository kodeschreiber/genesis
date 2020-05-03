import time

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

lol = Grid(2,5)
for pt in lol:
  lol.set(pt, int((time.time() - int(time.time())) * 1000000))

print(lol.__dict__)
