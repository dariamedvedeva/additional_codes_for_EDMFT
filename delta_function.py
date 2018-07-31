import numpy as np
import sys

# Structure
# 0 - name of a script
# 1 - beta
# 2 - (n) number of fermionic orbitals
# 3 .. 3 + n - Epsk
# .. end - Vk

Epsk = []
Vk = []

if len (sys.argv) > 1:
    beta = float(sys.argv[1])
    print "beta = ", beta
    print "++++++++++++++++++++++"
    n = int(sys.argv[2]) # n - number of orbitals
    for i in range (n):
        j = 4 + i
        Epsk.append(float(sys.argv[j - 1]))
        Vk.append(float(sys.argv[j + (n - 1)]))
    print Epsk
    print Vk
else:
    print ("Check the number of parameters: 1 - name of your script\n 2 - number of fermionic orbitals")
    print ("3 - ")
    print ("4 .. 4 + n - Epsk (energy of a level)\n from (4 + n) .. to the end - Vk (hybridization parameter) ")
    os.abort()

print "++++++++++++++++++++++"

f = open("delta.dat", "w")

cmplx_one = complex (0.0, 1.0)

for i in range(20):
    omega = (2 * i + 1) * cmplx_one * np.pi / beta
    value = 0.0
    for num in range(n):
        value += Vk[num]**2 /(omega - Epsk[num])
    print omega.imag, value.imag
f.close()
