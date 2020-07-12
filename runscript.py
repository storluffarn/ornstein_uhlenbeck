from os import system

#spring = 1.0e10
#t0 = 5e-11
#
##arr = []
##for k in range (-15,15,5):
##    arr.append(10**k)
#
#arr = [1e-10]
#
##arr = [1]
#
#for k in arr :
#    system("g++ -std=c++17 ornuhlimage.cpp -o ornuhlimage -O3 && ./ornuhlimage {} {}".format(k*spring,t0))
#    system("python plots.py")
    
#system("g++ -std=c++17 ornuhlimage.cpp -o ornuhlimage -O3 && ./ornuhlimage")
system("g++ -std=c++17 -fopenmp ornuhlimage.cpp -o ornuhlimage -O3 && ./ornuhlimage")
#system("python plots.py")
