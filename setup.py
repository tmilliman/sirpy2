# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='sirpy2',
    version='0.1.0',
    description='Package to read .sir files with python',
    long_description=readme,
    author='Thomas Milliman',
    author_email='thomas.milliman@unh.edu',
    url='https://github.com/tmilliman/sirpy2',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)

