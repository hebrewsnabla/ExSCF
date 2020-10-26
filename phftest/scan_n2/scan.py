import matplotlib.pyplot as plt
import sys

#data = sys.argv[1]
with open('rhf_e' + '.dat', 'r') as f:
    curvedata = f.readlines()
curve = []
for i in curvedata:
   curve.append(float(i.strip('\n')))

with open('uhf_e' + '.dat', 'r') as f:
    curvedata_u = f.readlines()
curve_u = []
for i in curvedata_u:
   curve_u.append(float(i.strip('\n')))

with open('suhf_e' + '.dat', 'r') as f:
    curvedata_su = f.readlines()
curve_su = []
for i in curvedata_su:
   curve_su.append(float(i.strip('\n')))

r = [#0.5, 0.6, 0.65, 0.7, 
    0.75, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.7, 2.0, 2.2, 2.5]
ru = [0.75, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]


plt.figure()
rhf = plt.plot(r, curve, label='rhf')
uhf = plt.plot(ru, curve_u, label='uhf')
suhf = plt.plot(r, curve_su, label='suhf')
plt.legend()
plt.savefig('hf.png')
