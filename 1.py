import random
from operator import itemgetter
import matplotlib.pyplot as plt
import math

CXPB = 0.8 # 交叉概率
MUPB = 0.2 # 变异概率
NGEN = 1000 # 迭代次数
popsize = 100 # 种群规模

d = 30 # 变量维度
encoding_size = 10 # 编码长度，变量的某一维用10位二进制表示
'''
bi 01串的数组
返回映射到0-1之间的浮点数
'''
def binary_to_float(bi):
    str1 = "".join(bi) # 先将List变成str
    number = int(str1, 2) # 根据str转为整数
    number /= 1024
    return number
'''
bi 一个个体的基因, 长度为编码长度乘以维度
返回对应的浮点数数组,长度为维度
'''
def decode_(bi):
    res = []
    length = len(bi)
    size = math.floor(length/d)
    temp = [bi[i:i+size] for i in range(0, length, size)]
    for i in range(len(temp)):
        res.append(binary_to_float(temp[i]))
    return res


class NSGA():
    def __init__(self):
        '''
        先初始化种群
        pop 种群
        Q 交叉变异后的种群
        R pop和Q的集合，父代种群和子代种群的集合
        每个个体以字典形式存储 gene表示二进制编码,
        '''
        pop = [] # 种群
        for i in range(popsize):
            individual = [] # 存储个体的01编码
            for j in range(encoding_size*d):
                individual.append(random.choice('01'))
            f1 = self.func1(decode_(individual)) # 计算个体的函数值
            f2 = self.func2(decode_(individual))
            pop.append({'gene':individual, 'f1':f1, 'f2':f2})
        self.pop = pop
        self.Q = []
        self.R = []

    def func1(self, x): # 计算目标函数1
        res = x[0]
        return res
        # return 1./res

    def func2(self, x): # 目标函数2
        sum = 0.0
        for i in range(d-1):
            sum += x[i+1]
        sum *= 9
        g = (1 + sum/(d-1))
        res = g * (1 - math.sqrt(x[0]/g))
        return res
        # return 1./res

    '''
    交叉
    off1, off2 两个个体
    '''
    def cross(self, off1, off2):
        gen1 = off1['gene']
        gen2 = off2['gene']

        pos1 = random.randrange(1, encoding_size*d)
        pos2 = random.randrange(1, encoding_size*d)

        newOff1 = []
        newOff2 = []
        '''
        双点交叉
        '''
        for i in range(encoding_size*d):
            if min(pos1, pos2) <= i < max(pos1, pos2):
                newOff1.append(gen1[i])
                newOff2.append(gen2[i])
            else:
                newOff1.append(gen2[i])
                newOff2.append(gen1[i])
        return {
            'gene': newOff1,
            'f1' : self.func1(decode_(newOff1)),
            'f2' : self.func2(decode_(newOff2))
        },{
            'gene' : newOff2,
            'f1' : self.func1(decode_(newOff2)),
            'f2': self.func2(decode_(newOff2))
        }

    '''
    变异
    '''
    def mut(self, individual):
        # 每个维度都变异
        idx = [0 for _ in range(d)]
        for i in range(d):
            idx[i] = random.randrange(i*encoding_size, (i+1)*encoding_size)
        for i in range(d):
            individual['gene'][idx[i]] = random.choice('01')
        individual['f1'] = self.func1(decode_(individual['gene']))
        individual['f2'] = self.func2(decode_(individual['gene']))
        return individual

    '''
    定义支配关系，作为排序的依据
    p的函数值都不大于q，并且至少有一维小于q，则p支配q
    注意这里的写法与目标函数的个数有关
    '''
    def dominate(self, p, q): # p支配q
        pv1 = self.R[p]['f1']
        pv2 = self.R[p]['f2']

        qv1 = self.R[q]['f1']
        qv2 = self.R[q]['f2']
        if (pv1 <= qv1) and (pv2 <= qv2):
            if (pv1 == qv1) and (pv2 == qv2):
                return False
            else:
                return True
        else:
            return False
    '''
    快速非支配集排序
    '''
    def fast_non_dominate_sort(self, R):
        Sp = [[] for i in range(len(R))] # 被个体p支配的个体集合
        np = [0 for i in range(len(R))] # 支配个体p的个体数量
        rank = [0 for i in range(len(R))] # 每个个体的rank
        F = [[]] # 层级

        for p in range(len(R)):
            Sp[p] = []
            np[p] = 0
            for q in range(len(R)):
                if self.dominate(p, q):
                    Sp[p].append(q)
                elif self.dominate(q, p):
                    np[p] += 1
            if np[p] ==0:
                rank[p] = 0
                F[0].append(p)
        i = 0
        while(len(F[i]) != 0):
            Q = []
            for p in F[i]:
                for q in Sp[p]:
                    np[q] = np[q] -1
                    if np[q] == 0:
                        rank[q] = i+1
                        Q.append(q)
            i = i+1
            if len(Q) != 0:
                F.append(Q)
            else:
                break
        return F
    '''
    计算群体L的拥挤距离
    '''
    def crowd(self, L):
        l = len(L)
        dis = [0 for _ in range(l)]
        for i in range(2):
            f = 'f'+str(i+1)
            M = sorted(L, key=itemgetter(f), reverse=True) # 从大到小排序
            dis[0] = dis[l-1] = 100000 # 两边的距离为无穷大
            for i in range(1, l-1):
                if M[0][f] - M[l-1][f] != 0:
                    dis[i] += (M[i+1][f] - M[i-1][f])/(M[0][f] - M[l-1][f])
        return dis



    def main(self):
        for g in range(NGEN):
            Q = []
            for i in range(int(popsize/2)):
                off1 = self.pop[i]
                off2 = self.pop[i+1]
                if random.random() <= CXPB:
                    off1, off2 = self.cross(off1, off2)
                if random.random() <= MUPB:
                    off1 = self.mut(off1)
                    off2 = self.mut(off2)
                Q.append(off1)
                Q.append(off2)
            self.Q = Q

            # P, Q加入到一起成为R
            self.R = []
            for i in range(popsize):
                self.R.append(self.pop[i])
                self.R.append(self.Q[i])
            FRONT = self.fast_non_dominate_sort(self.R)

            # 需要从R中筛选出P群体
            nextP = []
            for i in range(len(FRONT)):
                if len(nextP) + len(FRONT[i]) <= popsize:
                    for j in FRONT[i]:
                        nextP.append(self.R[j])
                    # print([self.R[j] for j in FRONT[i]])
                    # nextP.extend(self.R[j] for j in FRONT[i])
                else:
                    # 还需要popsize-len(nextP)个
                    L = [self.R[j] for j in FRONT[i]]
                    index = [j for j in FRONT[i]]
                    dis = self.crowd(L)
                    temp = []
                    for j in range(len(dis)):
                        temp.append({
                            'index' : index[j],
                            'crowd' : dis[j]
                        })
                    res = sorted(temp, key=itemgetter('crowd'), reverse=True) # 拥挤距离从小到大排序
                    rest = popsize-len(nextP)
                    for k in range(rest):
                        nextP.append(self.R[res[k]['index']])
                    break
            self.pop = nextP

            # if g == NGEN-1:
            print('iter:', g)
            res = [self.R[j] for j in FRONT[0]]
            print("x:{}, y1:{}, y2:{}".format(decode_(res[0]['gene']), res[0]['f1'], res[0]['f2'])  )

if __name__ == '__main__':
    nsga = NSGA()
    nsga.main()
