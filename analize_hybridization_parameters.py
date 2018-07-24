import numpy as np

if len (sys.argv) > 1:
    filename = sys.argv[1]
else:
    print ("Enter filename")
    os.abort()

f = open(filename, "r")
f_lambda = open("lambda_for_iterations.dat", "w")
lines = f.readlines()
count = 0
array = []
for line in lines:
    if count < 3:
        array.append(float(line.split()[0]))
    if count == 2:
        f_lambda.write(str(array[0]))
    	f_lambda.write(' ')
	f_lambda.write(str(array[1]))
	f_lambda.write(' ')
	f_lambda.write(str(array[2]))
	f_lambda.write('\n')
    count += 1

    if count >= 6:
        array[:] = []
        count = 0
    else:
        continue

f.close()
f_lambda.close()
