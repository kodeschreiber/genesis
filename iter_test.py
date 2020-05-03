import time

class Grid:
  def __init__(self, width, length):
    self.__xpos = 0
    self.__ypos = 0
    self.__fully_init = False
    self.__is_iter = False
    self.__width = width
    self.__length = length
    self.__grid_data = list()
    for i in range(width):
      self.__grid_data.append(list())
      for j in range(length):
        self.__grid_data[i].append(None)

  def __iter__(self):
    self.__xpos = 0
    self.__ypos = 0
    self.__is_iter = True
    return self

  def __next__(self):
    res = None
    if not self.__fully_init:
      res = (self.__xpos, self.__ypos)
    else:
      res = self.__grid_data[self.__xpos][self.__ypos]

    if self.__xpos < self.__width - 1:
      self.__xpos += 1
    else:
      self.__xpos = 0
      self.__ypos += 1
      if self.__ypos >= self.__length:
        if not self.__fully_init:
          self.__fully_init = True
        self.__is_iter = False
        raise StopIteration
    return res

  def __getitem__(self, idx):
    return self.__grid_data[idx]

  def __setitem__(self, idx, val):
    self.__grid_data[idx] = val

lol = Grid(3,2)
for x,y in lol:
  lol[x][y] = time.time()

print(lol.__dict__)
