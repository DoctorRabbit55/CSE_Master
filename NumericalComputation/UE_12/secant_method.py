import numpy as np
import matplotlib.pyplot as plt

def newton(x, f, df):
    return x - f(x)/df(x)

def f(x):
    return 2 - x**2 - np.exp(x)

def df(x):
    return -2*x - np.exp(x)

def secant(x1, x0, f, df):
    f1 = f(x1)
    f0 = f(x0)

    if f1-f0 == 0:
        f_diff = df(x0)
    else:
        f_diff = f1-f0
    
    return x1 - (x1-x0)/(f_diff)*f(x1)


N = 8
x_newton = 2.5
x_secant = 2.5
x_gt = 0.53727444917385660426

x_list = []
x_list.append(x_secant)
x_list.append(newton(x_secant, f, df))

error_list_newton = []
error_list_secant = []

# comparison between secand and newton method

for n in range(N):

    x_newton = newton(x_newton, f, df)
    error_list_newton.append(abs(x_newton-x_gt))

    x_list.append(secant(x_list[-1], x_list[-2], f, df))   

    error_list_secant.append(abs(x_list[-1]-x_gt))

plt.semilogy(range(N), error_list_newton, label='newton')
plt.semilogy(range(N), error_list_secant, label='secant')
plt.xlabel('iterations')
plt.ylabel('absolute error')
plt.legend()
plt.show()

# convergence order
p_list = []
for n in range(2,N):
    p_list.append(np.log(abs(x_list[n+1]-x_gt)) / np.log(abs(x_list[n]-x_gt)))

print('p_n = ',p_list)
plt.plot(range(2,N), p_list)
plt.xlabel('iterations')
plt.ylabel('convergence order')
plt.show()

# comparison in regard to efficiency

def compute_accuracy(x, x_gt):

    accuracy = 0
    
    for i in range(0,16):
        
        if int(x) == int(x_gt):
            accuracy += 1
        else:
            break

        x = x * 10
        x_gt = x_gt * 10

    return accuracy

print(error_list_secant)

accuracy_newton = []
accuracy_secant = []

for i in range(len(error_list_newton)):
    accuracy_newton.append(compute_accuracy(error_list_newton[i], 0) / float((i+1)*3))
    #accuracy_newton.append(error_list_newton[i] / float((i+1)*3))  


for i in range(len(error_list_secant)):
    accuracy_secant.append(compute_accuracy(error_list_secant[i], 0) / float((i+1)*2))
    #accuracy_secant.append(error_list_secant[i] / float((i+1)*3))  

plt.plot(range(len(accuracy_newton)), accuracy_newton, label='newton')
plt.plot(range(len(accuracy_secant)), accuracy_secant, label='secant')
plt.xlabel('iteration')
plt.ylabel('efficiency')
plt.legend()
plt.show()
