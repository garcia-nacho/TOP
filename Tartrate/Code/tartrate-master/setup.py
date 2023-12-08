#from setuptools import setup
from distutils.core import setup
from distutils.util import convert_path

versionholder = {}
ver_path = convert_path('tartrate/__init__.py')
with open(ver_path) as init_file:
    exec(init_file.read(), versionholder)

setup(name='tartrate',
    version=versionholder['__version__'],
    description='Script to differentiate typhoid and non-typhoid Salmonella by analysis of the STM3356 ORF',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],
    keywords=['tartrate','Salmonella','typhoid','paratyphoid','STM3356'],
    author='Ola Brynildsrud',
    author_email='ola.brynildsrud@fhi.no',
    install_requires=[
        'biopython'],
    entry_points={
        'console_scripts': ['tartrate=tartrate.tartrate:main']
    },
    packages=['tartrate'],
    package_dir={'tartrate': 'tartrate'},
    package_data={'tartrate': ['tartrate/db/*']},
    include_package_data=True
)