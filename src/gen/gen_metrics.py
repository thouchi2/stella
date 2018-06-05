import numpy as np
from mako.template import Template
import itertools
from gen_util import *


class metric_gen:

    def __init__(self):
        names = get_names()
        self.iname = names['i']
        self.pname = names['p']
        self.cname = names['c']
        self.lname = names['l']
        self.rep = ''


    def gen_init(self):
        ret = ''

        ret += 	'enum {c=0,xh=1,yh=2,zh=3};\n'

        tpl_metdef = Template('PetscScalar ***vec__${p}_${l};\ndouble ${p}_${l}[4];\nvec__${p}_${l} = jac_v[met->t3map[${i}][${j}]];\n')
        for i in range(3):
            for j in range(3):
                ret += tpl_metdef.render(p=self.pname[i], l=self.lname[j], i=i, j=j)

        for i in range(3):
            for j in range(3):
                if i <= j:
                    ret += 'double g_%d%d;\n' % (i+1, j+1)
                    ret += 'double g%d%d;\n' % (i+1, j+1)

        ret += 'double g;\n'
        ret += 'double pt0[3], pt1[3], pt2[3];\n'

        for i in range(3):
            ret += '%sbeg = %ss; %send = %ss + %sm; ' %(self.iname[i], self.pname[i],
                                                        self.iname[i], self.pname[i], self.pname[i])
        ret = ret[:-1] + '\n'

        for i in range(3):
            ret += 'if (%ss == 0) %sbeg++;\n' % (self.pname[i], self.iname[i])
            ret += 'if (%ss + %sm == ng%s) %send--;\n' % (self.pname[i], self.pname[i],
                                                          self.pname[i], self.iname[i])
        return ret


    def gen_store_metrics(self):
        '''Generate code to store metrics'''
        ret = ''

        ret += 'jac[k][j][i] = sqrt(g);\n';
        tpl_cpy = Template('vec__${p}_${l}[k][j][i] = ${p}_${l}[c];\n')
        for i in range(3):
            for j in range(3):
                ret += tpl_cpy.render(p=self.pname[j], l=self.lname[i])

        return ret


    def gen_deriv_cross(self, ddir, side=0):
        '''Generates grid metrics that influence cross terms in stencil.

        The metrics needed for cross terms are centered at grid points.

        Args:
            ddir: direction of difference
            side: indicates whether sided differences are needed
                  -1 for left, 1 for right, 0 for centered difference
        '''
        ret = ''
        idx0 = [0, 0, 0]
        idx1 = [0, 0, 0]
        idx2 = [0, 0, 0]

        if side == 0:
            idx0[ddir] = -1
            idx1[ddir] = 1
            for i in range(3):
                ret += assign(deriv(i, ddir),
                              gen('si[-1]') * coors(idx0, i) +
                              gen('si[1]') * coors(idx1, i))
        elif side == -1:
            idx1[ddir] = 1
            idx2[ddir] = 2
            for i in range(3):
                ret += assign(deriv(i, ddir),
                              gen('sm[0]') * coors(idx0, i) +
                              gen('sm[1]') * coors(idx1, i) +
                              gen('sm[2]') * coors(idx2, i))

        elif side == 1:
            idx1[ddir] = -1
            idx2[ddir] = -2
            for i in range(3):
                ret += assign(deriv(i, ddir),
                              gen('sp[0]') * coors(idx0, i) +
                              gen('sp[-1]') * coors(idx1, i) +
                              gen('sp[-2]') * coors(idx2, i))
        else:
            pass

        return ret


    def gen_deriv_coord(self, cdir, ddir, side=0):
        '''Generates metric terms that influence stencil coefficients in coordinate directions.

        These metric terms are centered at have points embedded in cdir lines.

        Args:
            cdir: half point is along the cdir direction
            ddir: direction of difference
            side: indicates whether sided differences are needed
                  -1 for left, 1 for right, 0 for centered difference
        '''
        ret = ''

        if cdir == ddir:
            idx0 = [0, 0, 0]
            idx1 = [0, 0, 0]
            idx1[ddir] = -1
            for i in range(3):
                ret += assign(deriv(i, ddir, cdir),
                              coors(idx0, i) - coors(idx1, i))
        else:
            ret += self.gen_interp_diff(cdir, ddir, side)

        return ret


    def gen_interp_diff(self, cdir, ddir, side=0):
        '''Generates differences needed at half points

        Since finite differences needed at half points when
        cdir != ddir require values at half points, this function
        first interpolates the needed points and then performs the
        difference.

        Args:
            cdir: half point is along the cdir direction
            ddir: direction of difference
            side: indicates whether sided differences are needed
                  -1 for left, 1 for right, 0 for centered difference
        '''
        ret = ''

        idx00 = [0, 0, 0]
        idx01 = [0, 0, 0]
        idx10 = [0, 0, 0]
        idx11 = [0, 0, 0]
        idx20 = [0, 0, 0]
        idx21 = [0, 0, 0]

        if side == 0:
            for l in range(3):
                if cdir == l:
                    idx00[l] = -1
                    idx10[l] = -1
                if ddir == l:
                    idx00[l] = -1
                    idx01[l] = -1
                    idx10[l] = 1
                    idx11[l] = 1

            for i in range(3):
                ret += assign(gen('pt0[%d]' %(i)),
                              gen('.5') * coors(idx00, i) +
                              gen('.5') * coors(idx01, i))
                ret += assign(gen('pt1[%d]' %(i)),
                              gen('.5') * coors(idx10, i) +
                              gen('.5') * coors(idx11, i))

            for i in range(3):
                ret += assign(deriv(i, ddir, cdir),
                              '(pt1[%d] - pt0[%d]) * .5' %(i, i))
        else:
            wt = []
            if side < 0:
                wt_name = 'sm'
            else:
                wt_name = 'sp'
            for l in range(3):
                wt.append('%s[%d]' %(wt_name,-1*side*l))
                if cdir == l:
                    idx00[l] = -1
                    idx10[l] = -1
                    idx20[l] = -1
                if ddir == l:
                    idx10[l] = -1*side
                    idx11[l] = -1*side
                    idx20[l] = -2*side
                    idx21[l] = -2*side
            for i in range(3):
                ret += assign(gen('pt0[%d]' %(i)),
                              gen('.5') * coors(idx00, i) +
                              gen('.5') * coors(idx01, i))
                ret += assign(gen('pt1[%d]' %(i)),
                              gen('.5') * coors(idx10, i) +
                              gen('.5') * coors(idx11, i))
                ret += assign(gen('pt2[%d]' %(i)),
                              gen('.5') * coors(idx20, i) +
                              gen('.5') * coors(idx21, i))
            for i in range(3):
                ret += assign(deriv(i, ddir, cdir),
                              '%s * pt0[%d] + %s * pt1[%d] + %s * pt2[%d]' %(wt[0], i, wt[1], i, wt[2], i))

        return ret


    def gen_cov(self, line_name):
        '''Generates code for covariant metric tensor'''
        ret = ''
        tpl = Template('g_${i+1}${j+1} = x_${lname[i]}[${line_name}]*x_${lname[j]}[${line_name}] + y_${lname[i]}[${line_name}]*y_${lname[j]}[${line_name}] + z_${lname[i]}[${line_name}]*z_${lname[j]}[${line_name}];\n')
        def gg(i,j):
            return tpl.render(i=i,j=j,line_name=line_name, lname=self.lname)

        for i in range(3):
            for j in range(3):
                if i <= j:
                    ret += gg(i, j)

        return ret


    def gen_det(self):
        '''Generates code for computing determinant of covariant metric tensor'''
        ret = ''
        ret += 'g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\n'

        return ret


    def gen_cont(self, i, j):
        '''Generates code for computing the contravariant tensor'''
        ret = ''

        if i == 0 and j == 0:
            ret += 'g11 = FDIV(1.0,g) * (g_22*g_33 - g_23*g_23);\n'
        if i == 1 and j == 1:
            ret += 'g22 = FDIV(1.0,g) * (g_11*g_33 - g_13*g_13);\n'
        if i == 2 and j == 2:
            ret += 'g33 = FDIV(1.0,g) * (g_11*g_22 - g_12*g_12);\n'
        if i == 0 and j == 1:
            ret += 'g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\n'
        if i == 0 and j == 2:
            ret += 'g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\n'
        if i == 1 and j == 2:
            ret += 'g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\n'

        return ret


    def gen_interior_coord(self, cdir):
        '''Generate interior metrics that influence stencil coefficients in coordinate directions.'''

        lbody = ''
        for deriv in range(3):
            lbody += self.gen_deriv_coord(cdir, deriv) + '\n'

        lbody += self.gen_cov(self.pname[cdir]+'h') + '\n'
        lbody += self.gen_det()
        lbody += self.gen_cont(cdir, cdir)

        idx = [0, 0, 0]
        idx[cdir] = -1
        lbody += assign(coef(idx, [0,0,0]),
                        gen('sqrt(g)') * cont(cdir+1, cdir+1))

        lbody = lbody[:-1]
        ends = ['iend', 'jend', 'kend']
        ends[cdir] = '%ss + %sm' % (self.pname[cdir], self.pname[cdir])
        return (loop('k', 'kbeg', ends[2],
                     loop('j', 'jbeg', ends[1],
                          loop('i', 'ibeg', ends[0], lbody))))



    def gen_set_cross(self, bnd=[0,0,0]):
        '''Generate code for assigning cross coefficients

        Args:
            bnd: indicates which faces this code touches
        '''
        ret = ''
        aop = '+'

        for dirs in [[-1,-1,0],
                     [1,-1,0],
                     [-1,0,-1],
                     [1,0,-1],
                     [0,-1,-1],
                     [0,1,-1]]:
            idx_cont = []
            if sum(dirs) == 0:
                aop = '-'
            else:
                aop = '+'
            for i in range(3):
                if dirs[i] !=0:
                    idx_cont.append(i+1)

            idx0 = [0, 0, 0]
            idx1 = [0, 0, 0]
            for i in range(3):
                idx0[i] = -1*dirs[i]
                if dirs[i] != 0:
                    break

            flag = False
            for i in range(3):
                if not flag and dirs[i] != 0:
                    flag = True
                    continue
                if flag:
                    idx1[i] = -1*dirs[i]

            for idx in [idx0, idx1]:
                safe = True
                for i in range(3):
                    if idx[i] == bnd[i] and idx[i] != 0:
                        safe = False
                if safe:
                    ret += assign(coef(dirs, idx),
                                  gen('.25') * gen('sqrt(g)') *
                                  cont(*idx_cont), aop)

        return ret


    def gen_interior_cross(self):
        ret = ''
        for i in range(3):
            for j in range(3):
                if i < j:
                    ret += self.gen_cont(i, j)

        ret += self.gen_set_cross()

        return ret


    def gen_interior(self):
        '''Generates metrics on interior of grid'''
        ret = '#define INTERIOR_METRICS_3D {\n'
        ret += self.gen_init()
        lbody = ''

        for i in range(3):
            lbody += self.gen_deriv_cross(i)

        lbody += self.gen_cov('c') + '\n'
        lbody += self.gen_det()
        lbody += self.gen_store_metrics()
        lbody += self.gen_interior_cross()

        ret += (loop('k', 'kbeg', 'kend',
                     loop('j', 'jbeg', 'jend',
                          loop('i', 'ibeg', 'iend', lbody)))) + '\n\n'

        for i in range(3):
            ret += self.gen_interior_coord(i) + '\n\n'

        ret += '}'
        return ret


    def gen_face_cross(self, fdir):
        ret = ''
        lbody = ''

        for dim in range(3):
            if dim == (abs(fdir) - 1):
                lbody += self.gen_deriv_cross(dim, fdir/abs(fdir))
            else:
                lbody += self.gen_deriv_cross(dim)
        lbody += self.gen_cov('c') + '\n'
        lbody += self.gen_det()
        lbody += self.gen_store_metrics()
        for i in range(3):
            for j in range(3):
                if i < j:
                    lbody += self.gen_cont(i,j)

        bnd = [0,0,0]
        bnd[abs(fdir)-1] = fdir/abs(fdir)
        lbody += self.gen_set_cross(bnd)

        ps = []
        for i in range(3):
            if (abs(fdir) - 1) != i:
                ps.append(self.iname[i])

        ret += loop(ps[1], '%sbeg' %(ps[1]), '%send' %(ps[1]),
                    loop(ps[0], '%sbeg' %(ps[0]), '%send' %(ps[0]), lbody))

        return ret


    def gen_face_coord(self, fdir):
        '''Generate code for metric terms influencing coefficients in coordinate directions

        Args:
            fdir: specifies which face (e.g. -1 is south face)
        '''
        ret = ''

        # compute metrics for each coordinate direction except the face
        for cdir in range(3):
            lbody = ''
            if abs(fdir)-1 != cdir:
                for dim in range(3):
                    if dim == (abs(fdir) - 1):
                        lbody += self.gen_deriv_coord(cdir, dim, fdir/abs(fdir))
                    else:
                        lbody += self.gen_deriv_coord(cdir, dim)

                lbody += self.gen_cov(self.pname[cdir]+'h') + '\n'
                lbody += self.gen_det()
                lbody += self.gen_cont(cdir,cdir)
                idx = [0, 0, 0]
                idx[cdir] = -1
                lbody += assign(coef(idx, [0,0,0]),
                                gen('sqrt(g)') * cont(abs(cdir)+1, abs(cdir)+1))

                ins = []
                ps = []
                ends = []
                for i in range(3):
                    if (abs(fdir) - 1) != i:
                        ins.append(self.iname[i])
                        ps.append(self.pname[i])
                        if i == cdir:
                            ends.append('%ss + %sm' %(ps[-1],ps[-1]))
                        else:
                            ends.append('%send' % (ins[-1]))

                ret += loop(ins[1], '%sbeg' %(ins[1]), ends[1],
                            loop(ins[0], '%sbeg' %(ins[0]), ends[0], lbody))


        return ret


    def gen_face(self, fdir):
        ret = ''
        ibody = ''
        if fdir < 0:
            ibody += assign(self.iname[abs(fdir)-1], "0")
        elif fdir > 0:
            ibody += assign(self.iname[abs(fdir)-1], "ng%s - 1" %(self.pname[abs(fdir)-1]))
        ibody += self.gen_face_coord(fdir)
        ibody += self.gen_face_cross(fdir)
        tpl = [Template('${p}s == 0'), Template('${p}s + ${p}m == ng${p}')]
        cond = ''
        if fdir < 0:
            cond += tpl[0].render(p=self.pname[abs(fdir)-1])
        elif fdir > 0:
            cond += tpl[1].render(p=self.pname[abs(fdir)-1])

        ret += genif(cond, ibody) + '\n\n'

        return ret


    def gen_corner_cross(self, dirs):
        ret = ''

        for dim in range(3):
            ret += self.gen_deriv_cross(dim, np.sign(dirs[dim]))
        ret += self.gen_cov('c')
        ret += self.gen_det()
        ret += self.gen_store_metrics()
        for i in range(3):
            for j in range(3):
                if i < j:
                    ret += self.gen_cont(i,j)

        bnd = [np.sign(dirs[i]) for i in range(3)]
        ret += self.gen_set_cross(bnd)

        return ret


    def gen_corner(self, dirs):
        ret = ''
        ibody = ''

        tpl = [Template('${p}s == 0'), Template('${p}s + ${p}m == ng${p}')]
        conds = ['','', '']
        cnt = 0
        for i in range(3):
            if dirs[i] < 0:
                conds[cnt] += tpl[0].render(p=self.pname[i])
                cnt += 1
                ibody += assign(self.iname[i], "0")
            elif dirs[i] > 0:
                conds[cnt] += tpl[1].render(p=self.pname[i])
                cnt += 1
                ibody += assign(self.iname[i], 'ng%s - 1' %(self.pname[i]))
        cond = conds[0] + ' && ' + conds[1] + ' && ' + conds[2]

        ibody += self.gen_corner_cross(dirs)[:-1]
        ret += genif(cond, ibody)

        return ret


    def gen_edge_cross(self, dirs):
        ret = ''
        lbody = ''

        for dim in range(3):
            lbody += self.gen_deriv_cross(dim, np.sign(dirs[dim]))
        lbody += self.gen_cov('c')
        lbody += self.gen_det()
        lbody += self.gen_store_metrics()
        for i in range(3):
            for j in range(3):
                if i < j:
                    lbody += self.gen_cont(i,j)

        bnd = [np.sign(dirs[i]) for i in range(3)]
        lbody += self.gen_set_cross(bnd)

        for i in range(3):
            if dirs[i] == 0:
                ins = self.iname[i]

        ret += loop(ins, '%sbeg' %(ins), '%send' % (ins), lbody)

        return ret

    def gen_edge_coord(self, dirs):
        ret = ''
        lbody = ''

        for cdir in range(3):
            if dirs[cdir] == 0:
                for dim in range(3):
                    lbody += self.gen_deriv_coord(cdir, dim, np.sign(dirs[dim]))

                lbody += self.gen_cov(self.pname[cdir]+'h') + '\n'
                lbody += self.gen_det()
                lbody += self.gen_cont(cdir,cdir)
                idx = [0, 0, 0]
                idx[cdir] = -1
                lbody += assign(coef(idx, [0,0,0]),
                                gen('sqrt(g)') * cont(abs(cdir)+1, abs(cdir)+1))

                ret += loop(self.iname[cdir], '%sbeg' %(self.iname[cdir]), '%ss + %sm' %(self.pname[cdir],self.pname[cdir]),
                        lbody)

        return ret


    def gen_edge(self, dirs):
        ret = ''
        ibody = ''
        for i in range(3):
            if dirs[i] < 0:
                ibody += assign(self.iname[i], "0")
            elif dirs[i] > 0:
                ibody += assign(self.iname[i], "ng%s - 1" % (self.pname[i]))
        ibody += self.gen_edge_coord(dirs) + '\n\n'
        ibody += self.gen_edge_cross(dirs)
        tpl = [Template('${p}s == 0'), Template('${p}s + ${p}m == ng${p}')]
        conds = ['','']
        cnt = 0
        for i in range(3):
            if dirs[i] < 0:
                conds[cnt] += tpl[0].render(p=self.pname[i])
                cnt += 1
            elif dirs[i] > 0:
                conds[cnt] += tpl[1].render(p=self.pname[i])
                cnt += 1
        cond = conds[0] + ' && ' + conds[1]

        ret += genif(cond, ibody)

        return ret


    def gen_boundary(self):
        '''Generates metric terms on grid boundary.'''

        ret = '#define BOUNDARY_METRICS_3D {\n'
        ret += self.gen_init()
        for i in range(3):
            ret += self.gen_face(i+1)
            ret += self.gen_face(-(i+1))

        for dirs in itertools.product([1,-1], repeat=2):
            edges = []
            edges.append([0, dirs[0], dirs[1]])
            edges.append([dirs[0], 0, dirs[1]])
            edges.append([dirs[0], dirs[1], 0])
            for i in range(3):
                ret += self.gen_edge(edges[i]) + '\n\n'

        for dirs in itertools.product([-1,1], repeat=3):
            ret += self.gen_corner(dirs) + '\n\n'

        ret += '}'
        return ret


    def gen_div_macro(self):
        ret = '#include <float.h>\n'
        ret += '#define FDIV(num, den) (fabs(den) < DBL_MIN ? 0.0 : num / den)\n'
        return ret


    def generate_c(self):
        tmp = ''
        tmp += '#ifndef EF_GEN_H\n'
        tmp += '#define EF_GEN_H\n\n'
        tmp += self.gen_div_macro()
        for line in self.gen_interior().split('\n'):
            tmp += line.rstrip() + '\\\n'
        tmp += '\n\n'
        for line in self.gen_boundary().split('\n'):
            tmp += line.rstrip() + '\\\n'
        tmp += '\n'
        tmp += '#endif\n'

        self.generated_code = tmp

    def save_c(self, fname):
        f = open(fname, 'w')
        f.write(self.generated_code)
        f.close()

    def __str__(self):
        tmp = ''
        tmp += self.gen_div_macro()
        for line in self.gen_interior().split('\n'):
            tmp += line.rstrip() + '\n'
        tmp += '\n\n'
        for line in self.gen_boundary().split('\n'):
            tmp += line.rstrip() + '\n'
        tmp += '\n'

        return tmp


if __name__ == '__main__':
    tst = metric_gen()
    tst.generate_c()
    tst.save_c('../../../src/ef_gen.h')
