# %%
import concurrent.futures
import multiprocessing
import numpy as np
from scipy import integrate
import multiprocessing


def f(x):
    return np.log(1/(x**2))

def _process(max):
    return integrate.quad(f, 2, max)[0]

print('hello')

# %%
if __name__ == '__main__':
    max_list = np.arange(4, 40)
    thread_pool = multiprocessing.Pool(multiprocessing.cpu_count())
    res = thread_pool.map(_process, max_list)
    print(res)


# %%
