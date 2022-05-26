from setuptools import setup, find_packages, Extension

setup(
    name='mykmeanssp',
    version='0.1.0',
    author="-",
    install_requires=['invoke'],
    packages=find_packages(),
    description='',
    classifiers=[
        'Programming Language  :: Python ::  Implementation  ::  CPhython',
    ],
    ext_modules =[
        Extension(
            'mykmeanssp',
            ['kmeans.c'],
        ),
    ]
)





