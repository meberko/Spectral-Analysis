import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys,math

class SmartPick():
    def __init__(self):
        self.muno_sources = []
        self.potential_sources = []
        self.actual_sources = []
        with open('chandra_point_sources_table3.txt','r') as f:
            content = f.readlines()
            content = content[45:]
            content = [x.split() for x in content]
            for source in content:
                self.muno_sources.append(source)
        with open('all.srclist', 'r') as f:
            content = f.readlines()
            content = content[1:]
            content = [x.split() for x in content]
            for source in content:
                self.actual_sources.append(source[0])
                self.potential_sources.append(source[0][:8]+source[0][9:16])

    def find_sources(self):
        i=0
        with open('all_smarter.srclist', 'w') as f:
            f.write('; Master Source List\n')
            for src in self.muno_sources:
                if src[1] in self.potential_sources:
                    idx = self.potential_sources.index(src[1])
                    if float(src[12]) > 25:
                        print src[1]
                        f.write(self.actual_sources[idx]+'\n')

if __name__ == '__main__':
    print "****************************"
    print "*       Smart Pick         *"
    print "****************************"

    sp = SmartPick()
    sp.find_sources()
