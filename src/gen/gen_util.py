from mako.template import Template


def loop(i, beg, end, body):
    lbody = Template(filename='loop.txt')
    return lbody.render(ind=i, beg=beg, end=end, body=body)


def genif(cond, body):
    tpl = Template(filename='if.txt')
    return tpl.render(cond=cond, body=body)


def assign(lhs, rhs, op=''):
    return '%s %s= %s;\n' %(str(lhs), op, str(rhs))


def get_names():
    ret = {'i' : ['i', 'j', 'k'],
           'p' : ['x', 'y', 'z'],
           'l' : ['r', 's', 't'],
           'c' : ['W', 'S', 'B']}

    return ret


class gen(object):

    def __init__(self, rep=''):
        names = get_names()
        self.iname = names['i']
        self.pname = names['p']
        self.cname = names['c']
        self.lname = names['l']
        self.rep = rep


    def offset(self, offsets):
        for i in range(3):
            val = offsets[3-i-1]
            if val > 0:
                self.rep += '[%s+%d]' % (self.iname[3-i-1], val)
            elif val < 0:
                self.rep += '[%s%d]' % (self.iname[3-i-1], val)
            else:
                self.rep += '[%s]' % (self.iname[3-i-1])


    def __add__(self, other):
        return gen(str(self) + ' + ' + str(other))


    def __sub__(self, other):
        return gen(str(self) + ' - ' + str(other))


    def __mul__(self, other):
        return gen(str(self) + ' * ' + str(other))


    def __str__(self):
        return self.rep


class coors(gen):

    def __init__(self, offsets, coord):
        super(self.__class__, self).__init__()
        self.rep = 'coors'
        self.offset(offsets)
        self.rep += '.%s' % (self.pname[coord])


class coef(gen):

    def __init__(self, dirs, offsets):
        super(self.__class__, self).__init__()
        labels = ['', 'E', 'N', 'T', 'B', 'S', 'W']
        label = labels[dirs[0]] + labels[dirs[1]*2] + labels[dirs[2]*3]
        if dirs[0] != 0 and dirs[1] != 0:
            label = label[::-1]
        self.rep += 'coef[MET_%s]' % (label)
        self.offset(offsets)


class deriv(gen):

    def __init__(self, cdir, ddir, line=None):
        super(self.__class__, self).__init__()
        if line is None:
            self.rep += '%s_%s[c]' % (self.pname[cdir], self.lname[ddir])
        else:
            self.rep += '%s_%s[%sh]' % (self.pname[cdir], self.lname[ddir], self.pname[line])


class cov(gen):

    def __init__(self, i, j):
        super(self.__class__, self).__init__()
        self.rep += 'g_%d%d' %(i,j)


class cont(gen):

    def __init__(self, i, j):
        super(self.__class__, self).__init__()
        self.rep += 'g%d%d' %(i,j)
