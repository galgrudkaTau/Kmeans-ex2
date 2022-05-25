from setuptools import setup, find_packages, Extension

setup(
    name='mykmeanssp',
    version='0.1.0',
    author="-",
    install_requires=['invoke'],
    packages=find_packages(),
    description='',
    ext_modules =[
        Extension(
            'kmeanspp',
            ['kmeans.c'],
        ),
    ]
)



