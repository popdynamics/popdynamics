from setuptools import setup

setup(
    name='popdynamics',
    version='0.1',
    description='Framework for epidemiological modelling',
    url='http://github.com/popdynamics/popdynamics',
    author='Bosco Ho',
    author_email='apposite@gmail.com',
    license='MIT',
    py_modules=['basepop'],
    install_requires=[
      'numpy',
      'matplotlib',
      'scipy',
      'graphviz',
      'future'
    ],
    zip_safe=False)
