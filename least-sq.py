import numpy as np
from scipy.optimize import leastsq, curve_fit, least_squares
import matplotlib.pyplot as plt

np.random.seed(2021)

# データ生成
def makeData(a, b, c, num_points):
    rd = np.random.normal(loc=0, scale=0.5, size=num_points*10)
    x_test = np.linspace(0, 16, num_points * 10)
    y_test = a * x_test ** 2 + b * x_test + c + rd
    return x_test, y_test

def func(prm, x, y):
    a, b, c = prm[0], prm[1], prm[2]
    residual = y-(a*x**2+b*x+c)
    return residual

# データ生成
num_points = 16
x, y = makeData(0.1, -1.6, 11, num_points)
x0=np.array([0, 0, 0])

# 最小二乗計算 損失関数 linear
result = least_squares(func, x0, loss='linear', f_scale=0.1, args=(x, y))
a, b, c = result.x[0], result.x[1], result.x[2]
y_linear = a*x**2+b*x+c

# 最小二乗計算 損失関数 soft_l1 
result = least_squares(func, x0, loss='soft_l1', f_scale=0.1, args=(x, y))
a, b, c = result.x[0], result.x[1], result.x[2]
y_soft_l1 = a*x**2+b*x+c

# 最小二乗計算 損失関数 cauchy 
result = least_squares(func, x0, loss='cauchy', f_scale=0.1, args=(x, y))
a, b, c = result.x[0], result.x[1], result.x[2]
y_cauchy = a*x**2+b*x+c

# グラフ表示
plt.title('Test scipy.optimize.least_squares()')
plt.plot(x, y, 'bo', label='y-original')
plt.plot(x, y_linear, color='red', label='y_linear')
plt.plot(x, y_soft_l1, color='orange', label='y_soft_l1')
plt.plot(x, y_cauchy, color='green', label='y_cauchy')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc='best', fancybox=True, shadow=True)
plt.grid(True)
plt.savefig('graph.png', dpi=70)
plt.show()
