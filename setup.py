from setuptools import setup, find_packages, Extension

setup(
    name='mykmeanspp',
    version='0.1.0',
    author="Maya Aderka & Gal Grudka",
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

