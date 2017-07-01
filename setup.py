# encoding: utf-8

from setuptools import setup
from setuptools.command.install import install

setup(
    name='Blasperin',
    version='0.1',
    author='Franca AF (Alexander da Franca Fernandes)',
    author_email='alexander@francafernandes.com.br',
    license='BSD',
    description='Tool to BLAST AnEnDB Fasta protein files.',
    long_description='Tool to BLAST AnEnDB Fasta protein files.',
    scripts=['bin/blasperin'],
    packages=[ 'blasperin' ],
    platforms='Linux',
    url='http://bioinfoteam.fiocruz.br/blasperin'
)


