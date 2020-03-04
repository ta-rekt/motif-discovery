import timeit

cy = timeit.timeit('''onionDecompose_cython.run()''', setup='import onionDecompose_cython', number=10)
py = timeit.timeit('''onionDecompose_python.run()''', setup='import onionDecompose_python', number=10)

print(cy, py)
print('Cython is {}x faster'.format(py/cy))
